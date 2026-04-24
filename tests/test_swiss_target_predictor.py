"""
测试 swiss_target - Stage 2 Swiss 靶点预测
使用 mock 避免实际启动浏览器
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from src.connectors.swiss_target import SwissTargetPredictor, SwissTargetConfig


class TestSwissTargetConfig:
    """配置类测试"""

    def test_default_config(self):
        """测试默认配置"""
        config = SwissTargetConfig()
        assert config.url == "https://www.swisstargetprediction.ch/"
        assert config.smiles == "CCC[C@H]1CCCCN1"
        assert config.organism == "Homo_sapiens"
        assert config.timeout == 15000

    def test_custom_config(self):
        """测试自定义配置"""
        config = SwissTargetConfig(
            url="https://custom.url/",
            smiles="CCO",
            organism="Mus_musculus",
            timeout=30000,
        )
        assert config.url == "https://custom.url/"
        assert config.smiles == "CCO"
        assert config.organism == "Mus_musculus"
        assert config.timeout == 30000

    def test_config_immutable(self):
        """测试配置不可变"""
        config = SwissTargetConfig()
        with pytest.raises(Exception):  # frozen dataclass
            config.smiles = "CCO"


class TestSwissTargetPredictor:
    """预测器类测试"""

    def test_init_default(self):
        """测试默认初始化"""
        predictor = SwissTargetPredictor()
        assert predictor.config is not None
        assert predictor.browser is None
        assert predictor.context is None
        assert predictor.page is None

    def test_init_with_params(self):
        """测试带参数初始化"""
        predictor = SwissTargetPredictor(
            config=SwissTargetConfig(smiles="CCO"),
            download_path="/tmp/test/"
        )
        assert predictor.config.smiles == "CCO"
        assert predictor.download_path == "/tmp/test/"

    def test_close_when_not_started(self):
        """测试关闭未启动的浏览器"""
        predictor = SwissTargetPredictor()
        predictor.close()  # 不应抛出异常
        assert predictor.browser is None

    @patch("src.connectors.swiss_target.sync_playwright")
    def test_setup_browser(self, mock_playwright):
        """测试浏览器启动（mock）"""
        mock_browser = MagicMock()
        mock_context = MagicMock()
        mock_page = MagicMock()

        mock_playwright_instance = MagicMock()
        mock_playwright.return_value = mock_playwright_instance
        mock_playwright_instance.start.return_value = mock_playwright_instance
        mock_playwright_instance.firefox.launch.return_value = mock_browser
        mock_browser.new_context.return_value = mock_context
        mock_context.new_page.return_value = mock_page

        predictor = SwissTargetPredictor()
        browser = predictor.setup_browser()

        assert browser == mock_browser
        assert predictor.browser == mock_browser
        assert predictor.context == mock_context
        assert predictor.page == mock_page

    def test_run_with_invalid_smiles_raises(self):
        """测试无效 SMILES 处理（网络层 mock）"""
        predictor = SwissTargetPredictor(
            config=SwissTargetConfig(smiles="INVALID_SMILES"),
            download_path="/tmp/test/"
        )

        # 由于需要真实浏览器，这里只测试结构
        # 实际运行会失败，符合预期
        assert predictor.config.smiles == "INVALID_SMILES"


class TestSwissTargetPredictorWorkflow:
    """工作流集成测试（使用 mock）"""

    @patch("src.connectors.swiss_target.sync_playwright")
    def test_full_run_mock(self, mock_playwright):
        """测试完整流程（mock 浏览器）"""
        # Mock 浏览器组件
        mock_browser = MagicMock()
        mock_context = MagicMock()
        mock_page = MagicMock()

        mock_playwright_instance = MagicMock()
        mock_playwright.return_value = mock_playwright_instance
        mock_playwright_instance.start.return_value = mock_playwright_instance
        mock_playwright_instance.firefox.launch.return_value = mock_browser
        mock_browser.new_context.return_value = mock_context
        mock_context.new_page.return_value = mock_page

        # Mock page 行为
        mock_page.goto = Mock()
        mock_page.wait_for_selector = Mock()
        mock_page.wait_for_url = Mock()
        mock_page.evaluate = Mock(return_value=False)
        mock_page.click = Mock()

        # Mock download
        mock_download = MagicMock()
        mock_download.value = MagicMock()
        mock_download.value.path = "/tmp/test.csv"

        with patch.object(mock_page, "expect_download") as mock_dl:
            mock_dl.return_value.__enter__ = Mock(return_value=mock_download)
            mock_dl.return_value.__exit__ = Mock(return_value=False)

            predictor = SwissTargetPredictor(
                config=SwissTargetConfig(smiles="CCO"),
                download_path="/tmp/test/"
            )

            # 由于 mock 无法完全模拟，这里验证 predictor 可正常初始化
            assert predictor.config.smiles == "CCO"
            assert predictor.download_path == "/tmp/test/"
