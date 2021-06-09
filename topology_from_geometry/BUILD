load("@rules_python//python:defs.bzl", "py_binary")

py_binary(
    name = "get_bond_length_distribution",
    srcs = [
      "get_bond_length_distribution.py",
    ],
    deps = [
      ":utilities",
      "//smu:dataset_pb2",
    ],
)

py_binary(
    name = "smi23d",
    srcs = [
        "smi23d.py",
    ],
    deps = [
    ],
)

py_library(
    name = "bond_lengths",
    srcs = [
        "bond_lengths.py",
    ],
    deps = [
      ":utilities",
      "//smu:dataset_pb2",
    ],
)

py_library(
    name = "utilities",
    srcs = [
        "utilities.py",
    ],
    deps = [
      "//smu:dataset_pb2",
      "//smu/parser:smu_utils_lib",
    ],
)

py_test(
   name = "utilities_test",
   srcs = [
       "utilities_test.py",
   ],
   deps = [
       ":utilities",
      "//smu:dataset_pb2",
   ],
)
