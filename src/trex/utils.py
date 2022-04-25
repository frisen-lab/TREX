import logging


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super().format(record)
