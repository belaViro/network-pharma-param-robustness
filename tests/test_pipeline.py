"""
测试 scripts.run_pipeline - 流水线
使用 mock 避免实际网络请求和浏览器启动
"""

import pytest
import os
import tempfile
from unittest.mock import patch

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pipeline import stage1_fetch_all, stage2_predict_all
from scripts.run_pipeline import run_pipeline, run_from_smiles


class TestStage1FetchAll:
    """Stage 1 测试"""

    @patch("src.pipeline.stage1_fetch.fetch_smiles")
    def test_fetch_single_molecule(self, mock_fetch):
        """测试单个分子查询"""
        mock_fetch.return_value = {
            "name": "Aspirin",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "status": "成功"
        }

        results = stage1_fetch_all(["Aspirin"])
        assert len(results) == 1
        assert results[0]["name"] == "Aspirin"
        assert results[0]["status"] == "成功"

    @patch("src.pipeline.stage1_fetch.fetch_smiles")
    def test_fetch_multiple_molecules(self, mock_fetch):
        """测试多个分子查询"""
        mock_fetch.side_effect = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功"},
            {"name": "Ibuprofen", "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "status": "成功"},
            {"name": "Caffeine", "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C", "status": "成功"},
        ]

        results = stage1_fetch_all(["Aspirin", "Ibuprofen", "Caffeine"])
        assert len(results) == 3
        assert all(r["status"] == "成功" for r in results)

    @patch("src.pipeline.stage1_fetch.fetch_smiles")
    def test_fetch_with_failures(self, mock_fetch):
        """测试部分失败"""
        mock_fetch.side_effect = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功"},
            {"name": "Invalid", "smiles": None, "status": "未找到"},
        ]

        results = stage1_fetch_all(["Aspirin", "Invalid"])
        assert len(results) == 2
        assert results[0]["status"] == "成功"
        assert results[1]["status"] == "未找到"


class TestStage2PredictAll:
    """Stage 2 测试"""

    def _mock_predict_target(self, smiles, download_path, organism="Homo_sapiens", molecule_name=None):
        """mock 函数替代 predict_target（可 pickle）"""
        return {
            "smiles": smiles,
            "status": "成功",
            "csv": f"{download_path}/{smiles[:10]}_targets.csv",
            "name": molecule_name or smiles,
        }

    @patch("src.connectors.swiss_target.SwissTargetPredictor.run_batch")
    def test_predict_single_smiles(self, mock_run_batch):
        """测试单个 SMILES 预测"""
        mock_run_batch.return_value = [{
            "name": "Aspirin",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "status": "成功",
            "csv": "/tmp/Aspirin_targets.csv",
        }]
        stage1_results = [{
            "name": "Aspirin",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "status": "成功"
        }]
        results = stage2_predict_all(stage1_results, "/tmp", max_workers=1)
        assert len(results) == 1
        assert results[0]["status"] == "成功"
        mock_run_batch.assert_called_once()

    @patch("src.connectors.swiss_target.SwissTargetPredictor.run_batch")
    def test_predict_skips_failed_stage1(self, mock_run_batch):
        """测试跳过 Stage 1 失败的结果"""
        mock_run_batch.return_value = [{
            "name": "Aspirin",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "status": "成功",
            "csv": "/tmp/Aspirin_targets.csv",
        }]
        stage1_results = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功"},
            {"name": "Invalid", "smiles": None, "status": "未找到"},
        ]
        results = stage2_predict_all(stage1_results, "/tmp", max_workers=1)
        assert len(results) == 1
        called_smiles = mock_run_batch.call_args.kwargs["smiles_list"]
        assert called_smiles == [("CC(=O)OC1=CC=CC=C1C(=O)O", "Aspirin")]

    def test_predict_empty_list(self):
        """测试空列表"""
        results = stage2_predict_all([], "/tmp")
        assert len(results) == 0


class TestRunFromSmiles:
    """从 SMILES 直接运行测试"""

    @patch("scripts.run_pipeline.stage2_predict_all")
    def test_run_from_smiles_directly(self, mock_stage2):
        """测试从 SMILES 列表直接运行"""
        mock_stage2.return_value = [
            {"smiles": "CCO", "status": "成功", "csv": "/tmp/CCO_targets.csv"}
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            results = run_from_smiles(["CCO"], output_dir=tmpdir)
            assert len(results) == 1
            mock_stage2.assert_called_once()


class TestPipelineIntegration:
    """流水线集成测试（mock）"""

    @patch("scripts.run_pipeline.stage1_fetch_all")
    @patch("scripts.run_pipeline.stage2_predict_all")
    def test_run_pipeline_full(self, mock_stage2, mock_stage1):
        """测试完整流水线"""
        mock_stage1.return_value = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功"},
        ]
        mock_stage2.return_value = [
            {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功", "csv": "/tmp/Aspirin_targets.csv"},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_pipeline(["Aspirin"], output_dir=tmpdir)

            assert "stage1_results" in result
            assert "stage2_results" in result
            assert "report_path" in result
            assert os.path.exists(result["report_path"])

    @patch("scripts.run_pipeline.stage1_fetch_all")
    @patch("scripts.run_pipeline.stage2_predict_all")
    def test_pipeline_report_content(self, mock_stage2, mock_stage1):
        """测试报告内容"""
        mock_stage1.return_value = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功"},
        ]
        mock_stage2.return_value = [
            {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "status": "成功", "csv": "/tmp/Aspirin_targets.csv"},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_pipeline(["Aspirin"], output_dir=tmpdir)

            with open(result["report_path"], "r", encoding="utf-8") as f:
                content = f.read()
                assert "SwissTargetPrediction 靶点预测流水线报告" in content
                assert "Stage 1" in content
                assert "Stage 2" in content
                assert "Aspirin" in content
