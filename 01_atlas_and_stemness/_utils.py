import logging
import os
import sys
from datetime import datetime

import yaml


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def setup_logger(name, log_path):
    ensure_dir(os.path.dirname(log_path))
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers = []

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    logger.info("Log start: %s", datetime.now().isoformat())
    return logger


def read_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def require_path(path, logger, hint=None):
    if not os.path.exists(path):
        msg = f"Missing required path: {path}"
        if hint:
            msg += f" | {hint}"
        logger.error(msg)
        raise FileNotFoundError(msg)


def find_first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None
