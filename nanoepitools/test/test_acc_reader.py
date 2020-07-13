import sys
import nanoepitools.accessibility as acc

def main():
    filename = '/home/r933r/data/projects/nanopore/train_footprinting/pred/unmet/0/accessibility_prediction.tsv'
    acc.AccessibilityProfile.read(filename)

if __name__ == '__main__':
    main()
