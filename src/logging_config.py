"""
统一日志配置
所有入口脚本应首先调用 setup_logging()
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

from src.config import LOG_FORMAT, LOG_DATE_FORMAT, LOG_LEVEL, REPORTS_DIR


def setup_logging(
    level: str = None,
    log_to_file: bool = True,
    log_dir: Path = None,
) -> logging.Logger:
    """配置统一日志系统

    Args:
        level: 日志级别 (DEBUG/INFO/WARNING/ERROR)，默认从配置读取
        log_to_file: 是否同时写入日志文件
        log_dir: 日志文件目录

    Returns:
        根日志记录器
    """
    if level is None:
        level = LOG_LEVEL
    level = getattr(logging, level.upper(), logging.INFO)

    # 清除已有 handlers（避免重复）
    root = logging.getLogger()
    if root.handlers:
        root.handlers.clear()

    # 控制台 Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)
    console_handler.setFormatter(console_formatter)
    root.addHandler(console_handler)

    # 文件 Handler
    if log_to_file:
        if log_dir is None:
            log_dir = REPORTS_DIR
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"pipeline_{timestamp}.log"
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)
        file_handler.setFormatter(file_formatter)
        root.addHandler(file_handler)
        root.info(f"日志文件: {log_file}")

    root.setLevel(level)
    return root
