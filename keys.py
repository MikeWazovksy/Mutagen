import subprocess

# Options
puzzle = 71
flips = 5
cpu = 26

# Path to the text file containing the keys (one per line)
key_file = 'keys.txt'

def run_program(key):
    # Prepare command with the provided key
    command = [
        "./mutagen",             # Path to your program
        "-p", str(puzzle),       # Parameter for puzzle number
        "-t", str(cpu),          # Parameter for the number of processors
        "-f", str(flips),        # Parameter for the number of flips
        "-k", str(key)           # Provided key
    ]
    
    # Run the program and wait for completion
    result = subprocess.run(command, capture_output=True, text=True)

    # Output the entire program result (stdout and stderr)
    print(f"Output for key {key}:\n{result.stdout}")
    if result.stderr:
        print(f"Error output for key {key}:\n{result.stderr}")

    # Check the program output for a successful solution
    if "Solution saved" in result.stdout:
        print(f"WoW ! Congratulations !")
        return True  # Solution found
    elif "No solution found" in result.stdout:
        print(f"Trying again...")
        return False  # Solution not found
    else:
        print(f"Unexpected output: {result.stdout}")
        return False  # In case of unexpected conclusions, we also try again.

def load_keys_from_file():
    with open(key_file, 'r') as file:
        keys = [line.strip() for line in file.readlines() if line.strip()]
    return keys

# Load keys from the file
keys = load_keys_from_file()

# Iterate through each key in the list
for key in keys:
    if run_program(key):
        break  # If a solution is found, stop the loop
