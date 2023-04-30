def pytest_addoption(parser):
    parser.addoption(
        "--data_dir",
        action="append",
        default=[],
        help="directory containing test data",
    )


def pytest_generate_tests(metafunc):
    if "data_dir" in metafunc.fixturenames:
        metafunc.parametrize("data_dir",
            metafunc.config.getoption("data_dir"))