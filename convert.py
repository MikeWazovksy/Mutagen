# Path to the source file with hex values
input_file = 'MATCH.txt'
# Path to the file where the results will be written
output_file = 'keys.txt'

def convert_hex_to_int(hex_value):
    try:
        # Convert hex to number
        return int(hex_value, 16)
    except ValueError:
        # If there was an error during conversion, return None
        print(f"Error while converting: {hex_value}")
        return None

def process_hex_file():
    # Открываем файл с hex значениями для чтения
    with open(input_file, 'r') as infile:
        hex_values = infile.readlines()

    # Open the file with hex values ​​for reading
    with open(output_file, 'w') as outfile:
        for hex_value in hex_values:
            # Remove spaces and newline characters
            hex_value = hex_value.strip()
            if hex_value:  # Check that the string is not empty
                # Convert hex to number
                number = convert_hex_to_int(hex_value)
                if number is not None:
                    # Write the result to a new file
                    outfile.write(f"{number}\n")
                else:
                    # If there is an error during conversion, write down the error message
                    outfile.write(f"Ошибка: {hex_value}\n")

# Start processing
process_hex_file()
