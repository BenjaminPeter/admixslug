def test_version(script_runner):
    ret = script_runner.run(["admixslug", "--version"])
    assert ret.success
