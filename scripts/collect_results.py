import os
import json


def read_json_in_folders(root_dir):
    """
    Reads JSON files inside every folder under the given root directory.

    Parameters:
        root_dir (str): Path to the root directory
    """
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.lower().endswith(".json"):
                file_path = os.path.join(dirpath, file)
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        data = json.load(f)
                        print(f"\n--- Contents of {file_path} ---")
                        print(json.dumps(data, indent=4))
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")


# Example usage:
if __name__ == "__main__":
    root_directory = "/path/to/your/directory"
    read_json_in_folders(root_directory)
