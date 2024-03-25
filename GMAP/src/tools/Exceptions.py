import sys
from traceback import TracebackException as TbEx


class GMAPException(Exception):

    __module__ = "builtins"

    def __init__(self, message="", error_code="UNKNOWN", exception=False):
        self.message = message
        self.error_code = error_code
        self.exception = exception
        sys.tracebacklimit = 0

    def __str__(self):
        if self.exception:
            traceprint = TbEx.from_exception(self.exception).format()
            printtrace = "".join(traceprint)
            return (
                f"{self.message}\n"
                f"{printtrace}\n"
                "More information can be found in the documentation "
                f"user pages using the following error code: {self.error_code}"
            )

        return (
            f"{self.message}\n"
            "More information can be found in the documentation "
            f"user pages using the following error code: {self.error_code}"
        )

