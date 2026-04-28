"""
测试 pubchem.worker - Stage 1 PubChem 查询
"""

import pytest
from src.workers.pubchem_worker import fetch_smiles


class TestFetchSmiles:
    """Stage 1 Worker 测试"""

    @pytest.mark.network
    def test_fetch_aspirin(self):
        """测试获取 Aspirin 的 SMILES（需要网络）"""
        result = fetch_smiles("Aspirin")
        assert result["name"] == "Aspirin"
        assert result["status"] == "成功"
        assert result["smiles"] is not None
        assert len(result["smiles"]) > 5  # SMILES 长度合理

    @pytest.mark.network
    def test_fetch_caffeine(self):
        """测试获取 Caffeine 的 SMILES"""
        result = fetch_smiles("Caffeine")
        assert result["name"] == "Caffeine"
        assert result["status"] == "成功"
        assert result["smiles"] is not None

    @pytest.mark.network
    def test_fetch_invalid_molecule(self):
        """测试无效分子名称"""
        result = fetch_smiles("ThisIsNotARealMoleculeName12345")
        assert result["name"] == "ThisIsNotARealMoleculeName12345"
        assert result["status"] != "成功"
        assert result["smiles"] is None

    @pytest.mark.network
    def test_fetch_smiles_structure(self):
        """验证返回的 SMILES 基本格式"""
        result = fetch_smiles("Ibuprofen")
        if result["status"] == "成功":
            smiles = result["smiles"]
            assert isinstance(smiles, str)
            assert len(smiles) > 0
            # SMILES 应该只包含常见原子符号
            import re
            assert re.match(r'^[A-Za-z0-9+\-\(\)\[\]=#@\/%]+$', smiles)
