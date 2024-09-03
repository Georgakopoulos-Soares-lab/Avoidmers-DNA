run_test:
	echo "Running unit tests. Hang on..."
	pytest -x -s -vv tests

clean:
	echo "Invoking magic broom"
	find . -type d -name "__pycache__" -exec rm -rf {} \;

