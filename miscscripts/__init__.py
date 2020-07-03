from pathlib import Path
from typing import *

def checkdir(path:Path)->Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path
def main():
	pass


if __name__ == "__main__":
	main()