"""
PubChem API Client Module

RESTful API for PubChem compound data:
https://pubchem.ncbi.nlm.nih.gov/rest/pug/

Supports: SMILES search, compound properties, bioactivity, structure retrieval.
"""

import logging
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Union

import requests

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PubChemConfig:
    """PubChem API configuration."""
    base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    timeout: int = 30
    max_retries: int = 3
    user_agent: str = "Mozilla/5.0 (compatible; ResearchBot/1.0)"


class PubChemError(Exception):
    """PubChem API error."""
    pass


class PubChemClient:
    """PubChem REST API client."""

    def __init__(self, config: Optional[PubChemConfig] = None):
        self.config = config or PubChemConfig()
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": self.config.user_agent})

    def _get(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Make GET request with retry."""
        url = f"{self.config.base_url}/{endpoint}"
        for attempt in range(1, self.config.max_retries + 1):
            try:
                response = self._session.get(url, params=params, timeout=self.config.timeout)
                response.raise_for_status()
                return response.json()
            except requests.RequestException as e:
                logger.warning(f"Request attempt {attempt}/{self.config.max_retries} failed: {e}")
                if attempt == self.config.max_retries:
                    raise PubChemError(f"Failed after {self.config.max_retries} retries: {e}") from e
        return {}

    def _post(self, endpoint: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """Make POST request."""
        url = f"{self.config.base_url}/{endpoint}"
        try:
            response = self._session.post(url, json=data, timeout=self.config.timeout)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            raise PubChemError(f"POST request failed: {e}") from e

    # -------------------------------------------------------------------------
    # Compound Search
    # -------------------------------------------------------------------------

    def find_compounds_by_smiles(self, smiles: str) -> List[int]:
        """Find compound CIDs by SMILES (2D similarity search).

        Args:
            smiles: SMILES string.

        Returns:
            List of CIDs.
        """
        # SMILES goes in the URL path, not as a query param
        result = self._get(f"compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON")

        try:
            cid_list = result.get("IdentifierList", {}).get("CID", [])
            if isinstance(cid_list, int):
                return [cid_list]
            return cid_list or []
        except (KeyError, TypeError) as e:
            logger.warning(f"Unexpected response structure: {e}")
        return []

    def find_compounds_by_name(self, name: str) -> List[int]:
        """Find compound CIDs by name.

        Args:
            name: Compound name.

        Returns:
            List of CIDs.
        """
        result = self._get(f"compound/name/{name}/cids/json")
        try:
            cid_list = result.get("IdentifierList", {}).get("CID", [])
            if isinstance(cid_list, int):
                return [cid_list]
            return cid_list or []
        except (KeyError, TypeError):
            return []

    # -------------------------------------------------------------------------
    # Compound Properties
    # -------------------------------------------------------------------------

    def get_compound_properties(self, cid: Union[int, str]) -> Optional[Dict[str, Any]]:
        """Get compound properties (MW, formula, etc.).

        Args:
            cid: PubChem Compound ID.

        Returns:
            Dict with properties or None.
        """
        props = "XLogP,ExactMass,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,ConnectivitySMILES"
        result = self._get(f"compound/cid/{cid}/property/{props}/JSON")
        try:
            return result.get("PropertyTable", {}).get("Properties", [{}])[0]
        except (KeyError, IndexError):
            return None

    def get_compound_synonyms(self, cid: Union[int, str]) -> List[str]:
        """Get compound synonyms/aliases.

        Args:
            cid: PubChem Compound ID.

        Returns:
            List of synonyms.
        """
        result = self._get(f"compound/cid/{cid}/synonyms/JSON")
        try:
            return result.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
        except (KeyError, IndexError):
            return []

    def get_compound_names(self, cid: Union[int, str]) -> Dict[str, List[str]]:
        """Get all names for a compound.

        Args:
            cid: PubChem Compound ID.

        Returns:
            Dict with 'IUPACName', 'CanonicalSMILES', 'Synonym', etc.
        """
        result = self._get(f"compound/cid/{cid}/JSON")
        try:
            props = result.get("PropertyTable", {}).get("Properties", [{}])[0]
            synonyms = self.get_compound_synonyms(cid)
            return {
                "iupac_name": props.get("IUPACName", ""),
                "canonical_smiles": props.get("CanonicalSMILES", ""),
                "isomeric_smiles": props.get("IsomericSMILES", ""),
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": props.get("MolecularWeight", 0),
                "synonyms": synonyms,
            }
        except (KeyError, IndexError):
            return {}

    # -------------------------------------------------------------------------
    # BioActivity
    # -------------------------------------------------------------------------

    def get_bioactivity_by_cid(self, cid: Union[int, str], assay_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """Get bioactivity data for a compound.

        Args:
            cid: PubChem Compound ID.
            assay_type: Filter by assay type (e.g., 'IC50', 'Ki', 'EC50').

        Returns:
            List of bioactivity records.
        """
        endpoint = f"compound/cid/{cid}/assaysummary/JSON"
        result = self._get(endpoint)
        try:
            assays = result.get("AssaySummary", {}).get("Assay", [])
            if assay_type:
                assays = [a for a in assays if a.get("ActivityType") == assay_type]
            return assays
        except (KeyError, TypeError):
            return []

    def get_bioactivity_by_target(self, gene_name: str, activity_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """Get bioactivity data by target gene name.

        Args:
            gene_name: Target gene name (e.g., 'BRCA1').
            activity_type: Filter by activity type.

        Returns:
            List of bioactivity records.
        """
        result = self._get(f"target/gene/{gene_name}/assays/JSON")
        try:
            assays = result.get("AssayList", []).get("Assay", [])
            if activity_type:
                assays = [a for a in assays if a.get("ActivityType") == activity_type]
            return assays
        except (KeyError, TypeError):
            return []

    # -------------------------------------------------------------------------
    # Structure Retrieval
    # -------------------------------------------------------------------------

    def get_2d_sdf(self, cid: Union[int, str]) -> Optional[str]:
        """Get 2D structure in SDF format.

        Args:
            cid: PubChem Compound ID.

        Returns:
            SDF string or None.
        """
        url = f"{self.config.base_url}/compound/cid/{cid}/SDF"
        try:
            response = self._session.get(url, timeout=self.config.timeout)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logger.warning(f"Failed to get SDF for CID {cid}: {e}")
            return None

    def get_3d_sdf(self, cid: Union[int, str]) -> Optional[str]:
        """Get 3D conformer in SDF format.

        Args:
            cid: PubChem Compound ID.

        Returns:
            SDF string or None.
        """
        url = f"{self.config.base_url}/compound/cid/{cid}/record/JSON?record_type=3d"
        try:
            response = self._session.get(url, timeout=self.config.timeout)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logger.warning(f"Failed to get 3D SDF for CID {cid}: {e}")
            return None

    def get_svg_image(self, cid: Union[int, str], width: int = 300, height: int = 200) -> Optional[bytes]:
        """Get compound 2D structure as SVG image.

        Args:
            cid: PubChem Compound ID.
            width: Image width.
            height: Image height.

        Returns:
            SVG bytes or None.
        """
        url = f"{self.config.base_url}/compound/cid/{cid}/PNG"
        params = {"image_size": f"{width}x{height}"}
        try:
            response = self._session.get(url, params=params, timeout=self.config.timeout)
            response.raise_for_status()
            return response.content
        except requests.RequestException as e:
            logger.warning(f"Failed to get image for CID {cid}: {e}")
            return None

    # -------------------------------------------------------------------------
    # Batch Operations
    # -------------------------------------------------------------------------

    def get_compounds_batch(self, cids: List[Union[int, str]]) -> List[Dict[str, Any]]:
        """Get properties for multiple compounds.

        Args:
            cids: List of PubChem CIDs.

        Returns:
            List of property dicts.
        """
        if not cids:
            return []

        cid_str = ",".join(str(c) for c in cids)
        props = "XLogP,ExactMass,MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName"
        result = self._get(f"compound/cid/{cid_str}/property/{props}/JSON")

        try:
            return result.get("PropertyTable", {}).get("Properties", [])
        except (KeyError, TypeError):
            return []

    def get_cid_by_inchikey(self, inchikey: str) -> Optional[int]:
        """Get CID from InChIKey.

        Args:
            inchikey: InChIKey string.

        Returns:
            CID or None.
        """
        result = self._get(f"inchikey/{inchikey}/cids/JSON")
        try:
            cid_list = result.get("IdentifierList", {}).get("CID", [])
            return cid_list[0] if cid_list else None
        except (KeyError, IndexError, TypeError):
            return None

    # -------------------------------------------------------------------------
    # Convenience Methods
    # -------------------------------------------------------------------------

    def get_compound_by_name(
        self,
        name: str,
        include_synonyms: bool = False,
    ) -> List[Dict[str, Any]]:
        """Get compound info by name in one call.

        Args:
            name: Compound name.
            include_synonyms: Include synonyms list.

        Returns:
            List of dicts with CID, SMILES, formula, MW, etc.
        """
        cids = self.find_compounds_by_name(name)
        if not cids:
            return []

        results = []
        for cid in cids:
            props = self.get_compound_properties(cid)
            if not props:
                continue
            smiles = (
                props.get("IsomericSMILES")
                or props.get("CanonicalSMILES")
                or props.get("ConnectivitySMILES")
                or ""
            )
            mw = props.get("MolecularWeight")
            try:
                mw = float(mw) if mw else 0.0
            except (ValueError, TypeError):
                mw = 0.0
            info = {
                "cid": cid,
                "smiles": smiles,
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": mw,
            }
            if include_synonyms:
                info["synonyms"] = self.get_compound_synonyms(cid)
            results.append(info)
        return results

    def get_compound_by_smiles(
        self,
        smiles: str,
        include_synonyms: bool = False,
    ) -> List[Dict[str, Any]]:
        """Get compound info by SMILES in one call.

        Args:
            smiles: SMILES string.
            include_synonyms: Include synonyms list.

        Returns:
            List of dicts with CID, SMILES, formula, MW, etc.
        """
        cids = self.find_compounds_by_smiles(smiles)
        if not cids:
            return []

        results = []
        for cid in cids:
            props = self.get_compound_properties(cid)
            if not props:
                continue
            smiles_out = (
                props.get("IsomericSMILES")
                or props.get("CanonicalSMILES")
                or props.get("ConnectivitySMILES")
                or ""
            )
            mw = props.get("MolecularWeight")
            try:
                mw = float(mw) if mw else 0.0
            except (ValueError, TypeError):
                mw = 0.0
            info = {
                "cid": cid,
                "smiles": smiles_out,
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": mw,
            }
            if include_synonyms:
                info["synonyms"] = self.get_compound_synonyms(cid)
            results.append(info)
        return results

    # -------------------------------------------------------------------------
    # Utility
    # -------------------------------------------------------------------------

    def convert_smiles_to_cid(self, smiles: str) -> Optional[int]:
        """Convert SMILES to CID via 2D similarity search.

        Args:
            smiles: SMILES string.

        Returns:
            CID or None.
        """
        result = self._get(f"compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON")
        try:
            cid_list = result.get("IdentifierList", {}).get("CID", [])
            if isinstance(cid_list, int):
                return cid_list
            return cid_list[0] if cid_list else None
        except (KeyError, TypeError):
            return None

    def get_cid_from_pubchem_pug(self, smiles: str) -> Optional[int]:
        """Convert SMILES to CID via identity search.

        Args:
            smiles: SMILES string.

        Returns:
            CID or None.
        """
        result = self._get(f"compound/fastidentity/smiles/{smiles}/cids/JSON")
        try:
            cid_list = result.get("IdentifierList", {}).get("CID", [])
            if isinstance(cid_list, int):
                return cid_list
            return cid_list[0] if cid_list else None
        except (KeyError, TypeError):
            return None


__all__ = [
    "PubChemClient",
    "PubChemConfig",
    "PubChemError",
]
