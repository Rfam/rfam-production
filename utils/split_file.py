import argparse


def split_file(num_lines, filename):
    """
    Split a large file into smaller files based on number of lines per file
    :param num_lines: number of lines to split by
    :param filename: large file name
    """
    lines_per_file = num_lines
    small_file = None
    with open(filename) as big_file:
        for line_num, line in enumerate(big_file):
            if line_num % lines_per_file == 0:
                if small_file:
                    small_file.close()
                small_filename = 'small_file_{}.txt'.format(line_num + lines_per_file)
                small_file = open(small_filename, "w")
            small_file.write(line)
        if small_file:
            small_file.close()


def parse_arguments():
    """
    Parse the command line arguments
    :return: Argparse parser object
    """
    parser = argparse.ArgumentParser(description='Split a given file into smaller files based on the number of line')

    req_args = parser.add_argument_group("required arguments")
    req_args.add_argument('-n', help='number of lines per smaller file', required=True)
    req_args.add_argument('-f', help='a file to split', type=str, required=True)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    split_file(args.n, args.f)
