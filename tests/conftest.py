def pytest_addoption(parser):
    parser.addoption(
        "--data_dir",
        action="append",
        default=[],
        help="directory containing test data",
    )


def pytest_generate_tests(metafunc):
    option_value = metafunc.config.getoption("data_dir")
    if "data_dir" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("data_dir", [option_value])
