import numpy as np

def compare_files(file1, file2):
    """
    Compare the contents of two files containing numerical data.

    Parameters:
    file1 (str): Path to the first file.
    file2 (str): Path to the second file.

    Returns:
    bool: True if the files are identical, False otherwise.
    """
    try:
        # Load the data from both files
        data1 = np.loadtxt(file1)
        data2 = np.loadtxt(file2)

        # Compare the data
        if np.array_equal(data1, data2):
            print("The files are identical.")
            return True
        else:
            print("The files are different.")
            # Find the indices where the arrays differ
            differences = np.where(data1 != data2)
            for index in zip(*differences):
                print(f"Difference at index {index}: file1={data1[index]}, file2={data2[index]}")
            return False
            return False
    except Exception as e:
        print(f"An error occurred while comparing files: {e}")
        return False

# Example usage
file1 = '/Users/curtisj/git_working_copies/zazzie/src/sassie/calculate/asaxs/asaxs_methods/poutput.txt'
file2 = '/Users/curtisj/git_working_copies/asaxs_code/src/scat_label_and_sum/matlab/output.txt'
file1 = '/Users/curtisj/git_working_copies/zazzie/src/sassie/calculate/asaxs/asaxs_methods/poutput2.txt'
file2 = '/Users/curtisj/git_working_copies/asaxs_code/src/scat_label_and_sum/matlab/output2.txt'
print('file1 = ', file1)
print('file2 = ', file2)
compare_files(file1, file2)
