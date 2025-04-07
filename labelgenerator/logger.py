import sys
import textwrap
import time

import loguru
from labelgenerator import __version__

SCRIPT_START = time.perf_counter()


logger = loguru.logger
logger.remove()
logger.add(sys.stderr, format="{message}")


def get_header():
    return textwrap.dedent(
        f"PyDLG version {__version__} released by The Exelixis Lab\n"
        f"Developed by: Julia Haag\n"
        f"Latest version: https://github.com/tschuelia/LabelGenerator\n"
        f"Questions/problems/suggestions? Please open an issue on GitHub.\n",
    )


def log_runtime_information(message, log_runtime=True):
    if log_runtime:
        seconds = time.perf_counter() - SCRIPT_START
        fmt_time = time.strftime("%H:%M:%S", time.gmtime(seconds))
        time_string = f"[{fmt_time}] "
    else:
        time_string = ""
    logger.info(f"{time_string}{message}")


def log_runtime(total_runtime: int, runtime_name: str):
    hours, remainder = divmod(total_runtime, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        logger.info(
            f"{runtime_name}: {int(hours):02d}:{int(minutes):02d}:{seconds:02d} hours ({round(total_runtime)} seconds)."
        )
    elif minutes > 0:
        logger.info(
            f"{runtime_name}: {int(minutes):02d}:{int(seconds):02d} minutes ({round(total_runtime)} seconds)."
        )
    else:
        logger.info(f"{runtime_name}: {seconds:.2f} seconds.")
