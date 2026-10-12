"""``Storage`` protocol and file-based storage of projects, experiments and results.

See SPEC §8.1.
"""

from chromsim.storage.base import Storage
from chromsim.storage.file_storage import FileStorage

__all__ = ["FileStorage", "Storage"]
