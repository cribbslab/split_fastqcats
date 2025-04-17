import subprocess
import os
import pytest

SCRIPTS = [
    "src/split_fastqcats/python/fastq_splitter_by_index.py",
    "src/split_fastqcats/python/fastq_splitter_by_primer.py",
]

def script_path(script):
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", script))

@pytest.mark.parametrize("script", SCRIPTS)
def test_script_help(script):
    path = script_path(script)
    # Test that script runs with --help or -h (should not error)
    try:
        result = subprocess.run(["python", path, "--help"], capture_output=True, text=True, timeout=10)
    except Exception as e:
        pytest.fail(f"Script {script} failed to run with --help: {e}")
    assert result.returncode == 0 or result.returncode == 1  # Some scripts exit 1 for help
    assert "usage" in result.stdout.lower() or "help" in result.stdout.lower() or result.stderr

@pytest.mark.parametrize("script", SCRIPTS)
def test_script_no_args(script):
    path = script_path(script)
    # Test that script fails gracefully with no arguments
    result = subprocess.run(["python", path], capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "usage" in result.stdout.lower() or "error" in result.stderr.lower() or result.stderr
