import os

if __name__ == "__main__":

    dirs = [
        x
        for x in os.listdir(os.getcwd())
        if os.path.isdir(os.path.join(os.getcwd(), x))
    ]

    fp = open("fixed_seeds_new.log", "w")

    for dir in dirs:
        if os.path.exists(os.path.join(dir, "SEED_old")) and os.path.exists(
            os.path.join(dir, "SEED")
        ):
            fp.write(dir + "\n")
    fp.close()
