from chemical_calculator import determine_organic_peroxide
import argparse

def parse_arguments():
    """
    Parse command line arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str,
                        help='The product list csv file to apply the chemical calculator on.')
    parser.add_argument('--save_file', type=str, default='results.csv',
                        help='The file name to save the results with.')
    args = parser.parse_args()
    input_file = args.input_file
    save_file = args.save_file

    return input_file, save_file

def main():
    input_file, save_file = parse_arguments()

    #batch process
    df_modified = determine_organic_peroxide(input_file)
    df_modified.to_csv(save_file)

if __name__ == '__main__':
    main()
