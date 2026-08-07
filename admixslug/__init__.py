from .interface_slug import run_sfs
from .interface_io import bam2, do_ref
from .interface_slug import profile as profile_slug
import importlib

__version__ = importlib.metadata.version(__name__)
