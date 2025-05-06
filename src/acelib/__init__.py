from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ace-elispot")
except PackageNotFoundError:
    __version__ = "unknown"