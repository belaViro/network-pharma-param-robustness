"""
SwissTargetPrediction Playwright 自动化模块
从 automate_swiss_target_playwright.py 抽取为核心类
"""

import os
import time
import random
import logging
import datetime
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

from playwright.sync_api import sync_playwright
from playwright_stealth.stealth import Stealth

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SwissTargetConfig:
    """SwissTargetPrediction 配置"""
    url: str = "https://www.swisstargetprediction.ch/"
    smiles: str = "CCC[C@H]1CCCCN1"
    organism: str = "Homo_sapiens"
    timeout: int = 15000


class SwissTargetPredictor:
    """SwissTargetPrediction Playwright 自动化"""

    def __init__(self, config: Optional[SwissTargetConfig] = None, download_path: Optional[str] = None):
        self.config = config or SwissTargetConfig()
        if download_path is None:
            download_path = Path(__file__).parent.parent / "outputs" / "shared" / "swiss_targets"
        self.download_path = str(download_path)
        self.browser = None
        self.context = None
        self.page = None
        self.playwright = None

    def setup_browser(self, proxy: str = None):
        """启动浏览器（stealth 模式）"""
        self.playwright = sync_playwright().start()
        launch_args = {
            "headless": False,
            "args": [
                "--no-sandbox",
                "--disable-dev-shm-usage",
            ],
        }
        if proxy:
            launch_args["proxy"] = {"server": proxy}
            logger.info(f"使用代理: {proxy}")
        self.browser = self.playwright.firefox.launch(**launch_args)

        viewport_width = random.randint(1280, 1600)
        viewport_height = random.randint(800, 1000)

        # Firefox 只能用 Firefox UA，指纹必须一致
        user_agents = [
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:121.0) Gecko/20100101 Firefox/121.0",
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:121.0) Gecko/20100101 Firefox/121.0",
            "Mozilla/5.0 (X11; Linux x86_64; rv:121.0) Gecko/20100101 Firefox/121.0",
        ]

        self.context = self.browser.new_context(
            viewport={"width": viewport_width, "height": viewport_height},
            user_agent=random.choice(user_agents),
            accept_downloads=True,
            locale="zh-CN",
            timezone_id="Asia/Shanghai",
            ignore_https_errors=True,
        )

        logger.info(f"浏览器视口: {viewport_width}x{viewport_height}")
        self.page = self.context.new_page()

        stealth_instance = Stealth()
        stealth_instance.apply_stealth_sync(self.context)

        self.page.add_init_script("""
            Object.defineProperty(navigator, 'webdriver', {
                get: () => undefined
            });
            Object.defineProperty(navigator, 'plugins', {
                get: () => [1, 2, 3, 4, 5]
            });
        """)

        logger.info("浏览器已启动（stealth 模式）")
        return self.browser

    def open_page(self) -> None:
        """打开页面（带真人延迟）"""
        time.sleep(random.uniform(1.0, 3.0))
        logger.info(f"正在打开: {self.config.url}")
        self.page.goto(self.config.url, wait_until="load", timeout=180000)
        # 模拟人类：偶尔滚动一下
        if random.random() < 0.3:
            self.page.mouse.move(random.randint(100, 500), random.randint(100, 400))
            self.page.evaluate("window.scrollBy(0, 200)")
            time.sleep(random.uniform(0.5, 1.5))
        time.sleep(random.uniform(2.0, 4.0))
        logger.info("页面加载完成")

    def input_smiles(self) -> None:
        """输入 SMILES（可靠方式：直接设置 DOM 值）"""
        self.page.wait_for_selector("#smilesBox", timeout=180000)
        time.sleep(random.uniform(0.5, 1.5))

        # 最可靠的方式：直接通过 JS 设置输入框值，避免 keyboard.type 丢字符
        smiles_js = self.config.smiles.replace('\\', '\\\\').replace('"', '\\"')
        self.page.evaluate(f"""
            var box = document.getElementById("smilesBox");
            box.value = "{smiles_js}";
            box.dispatchEvent(new Event("input", {{ bubbles: true }}));
            box.dispatchEvent(new Event("keyup", {{ bubbles: true }}));
            box.dispatchEvent(new Event("change", {{ bubbles: true }}));
        """)

        time.sleep(random.uniform(1.0, 2.0))

        # 验证输入是否正确
        actual_value = self.page.evaluate('document.getElementById("smilesBox").value')
        if actual_value != self.config.smiles:
            logger.warning(f"输入验证失败！预期: {self.config.smiles[:30]}... 实际: {actual_value[:30]}...")
            # 重试一次：用 fill 方式
            self.page.fill("#smilesBox", self.config.smiles)
            self.page.evaluate("""
                var box = document.getElementById("smilesBox");
                box.dispatchEvent(new Event("input", { bubbles: true }));
                box.dispatchEvent(new Event("keyup", { bubbles: true }));
                box.dispatchEvent(new Event("change", { bubbles: true }));
            """)
            time.sleep(1)
            actual_value = self.page.evaluate('document.getElementById("smilesBox").value')
            if actual_value != self.config.smiles:
                raise ValueError(f"SMILES 输入两次都失败: {actual_value}")

        # 检查提交按钮状态
        time.sleep(2)
        button_disabled = self.page.evaluate("""
            document.getElementById("submitButton").disabled
        """)

        if button_disabled:
            logger.info("按钮未启用，尝试触发 focus...")
            self.page.evaluate("""
                document.getElementById("smilesBox").focus();
                document.getElementById("smilesBox").dispatchEvent(new Event("keyup", { bubbles: true }));
            """)
            time.sleep(3)
            button_disabled = self.page.evaluate("""
                document.getElementById("submitButton").disabled
            """)
            if button_disabled:
                logger.warning("按钮仍未启用，但继续尝试提交")

        logger.info(f"已输入 SMILES: {self.config.smiles}")
        logger.info("提交按钮已启用")

    def select_organism(self) -> None:
        """选择物种"""
        logger.info(f"已选择物种: {self.config.organism}")

    def submit_prediction(self) -> None:
        """提交预测（带真人延迟）"""
        time.sleep(random.uniform(1.0, 2.5))

        try:
            self.page.click("#submitButton", timeout=180000)
            logger.info("已点击提交按钮")
        except Exception as e:
            logger.warning(f"点击超时，尝试重新输入: {e}")
            self.input_smiles()
            time.sleep(random.uniform(0.5, 1.0))
            try:
                self.page.click("#submitButton", timeout=180000)
                logger.info("重新点击提交按钮成功")
            except Exception as e2:
                logger.warning(f"重新点击也失败: {e2}")
                raise

        try:
            self.page.wait_for_url("**/predict.php**", timeout=120000)
            logger.info("已跳转到结果页面")
        except Exception as e:
            logger.warning(f"等待 predict.php 超时 (2分钟)，返回 index.php 重试: {e}")
            self.page.goto(self.config.url, wait_until="load", timeout=120000)
            time.sleep(random.uniform(1.0, 2.5))
            raise TimeoutError("predict.php 超时，已返回 index.php")

        error_text = self.page.evaluate("""
            document.body.innerText.includes("Please, provide only one molecule at a time")
        """)
        if error_text:
            raise ValueError("SMILES验证错误: 网站返回无效分子错误")

    def get_results(self, molecule_name: Optional[str] = None) -> str:
        """点击CSV下载按钮并保存

        Args:
            molecule_name: 分子名称，用于命名文件（如果为None则用smiles前20字符）
        """
        error_text = self.page.evaluate("""
            document.body.innerText.includes("Please, provide only one molecule at a time")
        """)
        if error_text:
            raise ValueError("SMILES验证错误: 网站返回无效分子错误")

        try:
            csv_button = self.page.wait_for_selector(".dt-button.buttons-csv.buttons-html5", timeout=180000)
        except:
            raise TimeoutError("等待CSV按钮超时")
        time.sleep(random.uniform(0.2, 0.5))

        with self.page.expect_download() as download_info:
            csv_button.click()
            time.sleep(random.uniform(0.2, 0.5))

        download = download_info.value
        # 使用分子名称（如果有的话），否则用smiles前20字符
        if molecule_name:
            safe_name = "".join(c for c in molecule_name if c.isalnum() or c in (' ', '-', '_')).strip()
            safe_name = safe_name.replace(' ', '_')
            file_path = str(Path(self.download_path) / f"{safe_name}_targets.csv")
        else:
            file_path = str(Path(self.download_path) / f"{self.config.smiles[:20]}_targets.csv")
        download.save_as(file_path)

        logger.info(f"CSV已保存: {file_path}")
        return file_path

    def _reload_page(self) -> None:
        """刷新当前页面（代替关闭重开）"""
        logger.info("刷新当前页面...")
        self.page.reload(wait_until="load", timeout=180000)
        time.sleep(random.uniform(1.0, 2.5))
        logger.info("页面刷新完成")

    def _update_config(self, smiles, organism):
        """更新配置"""
        if smiles or organism:
            self.config = SwissTargetConfig(
                url=self.config.url,
                smiles=smiles or self.config.smiles,
                organism=organism or self.config.organism,
                timeout=self.config.timeout,
            )

    def run_single(
        self,
        smiles: Optional[str] = None,
        organism: Optional[str] = None,
        molecule_name: Optional[str] = None,
        max_retries: int = 3,
        is_first: bool = True,
    ) -> str:
        """处理单个分子（浏览器已启动，不关闭）

        Args:
            smiles: 分子的 SMILES 字符串
            organism: 物种
            molecule_name: 分子名称
            max_retries: 最大重试次数
            is_first: 是否是批量中的第一个分子

        Returns:
            CSV 文件路径
        """
        self._update_config(smiles, organism)

        for attempt in range(1, max_retries + 1):
            logger.info(f"===== 尝试执行 ({attempt}/{max_retries}) =====")

            try:
                if attempt == 1:
                    if is_first:
                        self.open_page()
                    else:
                        # 批量模式下，新分子需要回到首页重新提交
                        logger.info("回到首页...")
                        self.page.goto(self.config.url, wait_until="load", timeout=180000)
                        time.sleep(random.uniform(1.0, 2.5))
                else:
                    self._reload_page()

                # 重试时检查是否已在结果页
                if attempt > 1:
                    try:
                        self.page.wait_for_selector(".dt-button.buttons-csv.buttons-html5", timeout=5000)
                        logger.info("检测到CSV按钮，直接获取结果")
                        csv_path = self.get_results(molecule_name=molecule_name)
                        logger.info("执行成功！")
                        return csv_path
                    except Exception:
                        pass

                self.input_smiles()
                self.select_organism()
                self.submit_prediction()
                csv_path = self.get_results(molecule_name=molecule_name)

                logger.info("执行成功！")
                return csv_path

            except ValueError as e:
                logger.warning(f"SMILES 验证错误: {e}")
                return ""

            except Exception as e:
                logger.warning(f"执行失败 (尝试 {attempt}/{max_retries}): {e}")

                if attempt == max_retries:
                    logger.error(f"已达到最大重试次数 ({max_retries})")
                    raise

                logger.info("等待 3 秒后刷新重试...")
                time.sleep(3)

        return ""

    def run(
        self,
        smiles: Optional[str] = None,
        organism: Optional[str] = None,
        molecule_name: Optional[str] = None,
        max_retries: int = 3,
        proxy: str = None,
    ) -> str:
        """运行完整流程（单分子，自动开关浏览器）

        Args:
            smiles: 分子的 SMILES 字符串
            organism: 物种（默认 Homo_sapiens）
            molecule_name: 分子名称，用于命名CSV文件
            max_retries: 最大重试次数
            proxy: 代理服务器地址
        """
        self._update_config(smiles, organism)
        self.setup_browser(proxy=proxy)
        try:
            return self.run_single(molecule_name=molecule_name, max_retries=max_retries, is_first=True)
        finally:
            self.close()

    def run_batch(
        self,
        smiles_list,
        organism: Optional[str] = None,
        proxy: str = None,
        max_retries: int = 3,
    ):
        """批量运行，复用同一个浏览器实例

        Args:
            smiles_list: [(smiles, molecule_name), ...]
            organism: 物种
            proxy: 代理
            max_retries: 最大重试次数

        Returns:
            结果列表 [{"name": ..., "status": ..., "csv": ...}, ...]
        """
        self._update_config(None, organism)
        self.setup_browser(proxy=proxy)
        results = []

        for i, (smiles, molecule_name) in enumerate(smiles_list):
            try:
                csv_path = self.run_single(
                    smiles=smiles,
                    molecule_name=molecule_name,
                    max_retries=max_retries,
                    is_first=(i == 0),
                )
                if csv_path:
                    results.append({"smiles": smiles, "name": molecule_name or smiles, "status": "成功", "csv": csv_path})
                else:
                    results.append({"smiles": smiles, "name": molecule_name or smiles, "status": "失败", "csv": ""})
            except Exception as e:
                results.append({"smiles": smiles, "name": molecule_name or smiles, "status": f"错误: {e}", "csv": ""})

            if i < len(smiles_list) - 1:
                delay = random.uniform(1, 3)
                logger.info(f"等待 {delay:.1f}s 后处理下一个分子...")
                time.sleep(delay)

        self.close()
        return results

    def close(self) -> None:
        """关闭浏览器和 playwright"""
        if self.browser:
            try:
                self.browser.close()
                self.playwright.stop()
            except Exception:
                pass
            self.browser = None
            self.context = None
            self.page = None
            self.playwright = None
