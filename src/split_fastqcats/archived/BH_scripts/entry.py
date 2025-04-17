import os
import sys
import glob
import importlib.util
import src

def main(argv=None):

    # Get the command-line arguments
    argv = sys.argv

    # Path to look for pipelines
    path = os.path.abspath(os.path.dirname(src.__file__))
    print(f"Searching for script in: {path}")

    # If no arguments or help is requested, display the help message
    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print(globals()["__doc__"])
        return

    # Attempt to find and load the corresponding pipeline
    else:
        print(globals()["__doc__"])

        # Searching for the pipeline script
        pipelines = glob.glob(os.path.join(path, "script_*.py"))
        pipeline_names = [os.path.basename(x)[len("script_"):-len(".py")] for x in pipelines]

        if not pipelines:
            print("Error: No scripts found.")
            sys.exit(1)

        print(f"Available scripts: {', '.join(pipeline_names)}\n")

        # Check if the pipeline argument is provided
        if len(argv) < 2:
            print("Error: No script specified.")
            sys.exit(1)

        pipeline = argv[1]
        matching_pipelines = [p for p in pipelines if os.path.basename(p) == f"script_{pipeline}.py"]

        if not matching_pipelines:
            print(f"Error: script '{pipeline}' not found.")
            sys.exit(1)

        pipeline_script = matching_pipelines[0]
        print(f"Loading script: {pipeline_script}")

        # Remove 'split_fastqcats' from sys.argv
        del sys.argv[0]

        # Load the pipeline module and execute its main function
        try:
            spec = importlib.util.spec_from_file_location(pipeline, pipeline_script)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            print(f"Successfully loaded module: {module.__name__}")
            module.main(sys.argv)
        except FileNotFoundError:
            print(f"Error: The pipeline script '{pipeline_script}' was not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading or executing the pipeline '{pipeline}': {e}")
            sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())
