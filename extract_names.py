import sys

if __name__ == "__main__":
    for l in sys.stdin:
        if "chromosome" in l:
            sys.stdout.write(l[1:])