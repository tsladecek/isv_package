from argparse import ArgumentParser
import pandas as pd
from isv import isv


if __name__ == "__main__":
    parser = ArgumentParser("Interpretation of Structural Copy Number Variants")

    parser.add_argument("-i", "--input", required=True, help="data to be predicted", type=str)
    parser.add_argument("-o", "--output", required=False, help="where the output tsv will be saved",
                        type=str, default="./isv_predictions.tsv")
    parser.add_argument("-p", "--proba", required=False, action="store_true", help="Return probabilities")
    parser.add_argument("-sv", "--shapvalues", required=False, help="Calculate SHAP Values", action="store_true")

    args = parser.parse_args()

    if args.input.endswith('.tsv') or args.input.endswith('.bed'):
        bed = pd.read_csv(args.input, sep='\t')
    elif args.input.endswith('.csv'):
        bed = pd.read_csv(args.input)
    else:
        exit("Unknown File extension. Use '.tsv', '.bed' or '.csv'")

    final = isv(cnvs=bed, proba=args.proba, shap=args.shapvalues)

    final.to_csv(args.output, sep='\t', index=False)
    print(f"Results saved to {args.output}")
