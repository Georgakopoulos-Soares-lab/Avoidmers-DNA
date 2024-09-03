install:
	echo "Installing dependencies...Please make sure you have created a new environment. This script will not automatically create a new one."
	pip install -r requirements.txt
	mamba install -c bioconda bedtools seqkit

run_test:
	echo "Running unit tests. Hang on..."
	pytest -x -s -vv tests

clean:
	echo "Invoking magic broom"
	find . -type d -name "__pycache__" -exec rm -rf {} \;

run:
	echo "Initializing snakemake pipeline..."
	rm -r -I avoidmers/
	bash submit_zimin_snake.sh
	echo "Running tests..."
	make run_test
	echo "ALL PASSED? SO FAR SO GOOD!"

