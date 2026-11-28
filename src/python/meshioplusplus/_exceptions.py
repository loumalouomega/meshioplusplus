class ReadError(Exception):
    pass


class WriteError(Exception):
    pass


class CorruptionError(Exception):
    pass


class CancelledError(Exception):
    pass
