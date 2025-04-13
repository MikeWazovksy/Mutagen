import random
import subprocess

# Options
puzzle = 69
flips = 5
cpu = 26
exclude = 4

# Constants for generating random number
LOWER = 295147905179352825856
UPPER = 590295810358705651711

def run_program(random_key):
    # Prepare command with random value for -k
    command = [
        "./mutagen",             # Path to your program
        "-p", str(puzzle),       # Parameter for puzzle number
        "-t", str(cpu),          # Parameter for the number of processors
        "-f", str(flips),        # Parameter for the number of flips
        "-x", str(exclude),      # Unchanged bits
        "-k", str(random_key)    # Generated random key
    ]
    
    # Run the program and wait for completion
    result = subprocess.run(command, capture_output=True, text=True)

    # Output the entire program result (stdout and stderr)
    print(f"{result.stdout}")
    if result.stderr:
        print(f"Error output for key {random_key}:\n{result.stderr}")

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

def generate_random_key():
    while True:
        # Generate a random key within LOWER and UPPER
        random_key = random.randint(LOWER, UPPER)
        
        # Convert the random key to hexadecimal
        hex_key = hex(random_key)[2:]  # Remove the '0x' prefix

        # Check if the first character is "1"
        if hex_key[0] == '1':
            return random_key  # Return the valid key
    
# Main loop to find the solution
while True:
    random_key = generate_random_key()
    if run_program(random_key):
        break  # If a solution is found, exit the loop
