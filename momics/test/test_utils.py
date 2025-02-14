import pytest
from momics.utils import *

"""
Testing utilities for resource monitoring, mainly memory usage and RAM load.
"""


# def test_memory_load(monkeypatch):
#     """
#     Tests the memory_load function.
#     """

#     class MockVirtualMemory:
#         def __init__(self, used, total):
#             self.used = used
#             self.total = total

#     def mock_virtual_memory():
#         return MockVirtualMemory(used=8 * (1024**3), total=16 * (1024**3))

#     monkeypatch.setattr("psutil.virtual_memory", mock_virtual_memory)

#     used_gb, total_gb = memory_load()

#     # Check if the used and total memory are correctly calculated
#     assert used_gb == 8, f"Expected used memory to be 8 GB, but got {used_gb} GB"
#     assert total_gb == 16, f"Expected total memory to be 16 GB, but got {total_gb} GB"


# def test_memory_usage(monkeypatch):
#     """
#     Tests the memory_usage function.
#     """

#     class MockVirtualMemory:
#         def __init__(self, percent):
#             self.percent = percent

#     def mock_virtual_memory():
#         return MockVirtualMemory(percent=50.0)

#     monkeypatch.setattr("psutil.virtual_memory", mock_virtual_memory)

#     percent_used = memory_usage()

#     # Check if the percentage of used memory is correctly calculated
#     assert (
#         percent_used == 50.0
#     ), f"Expected used memory percentage to be 50.0%, but got {percent_used}%"
