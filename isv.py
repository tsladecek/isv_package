from argparse import ArgumentParser
import pandas as pd
import sys

from isv import predict, shap_values
from isv.scripts.constants import HUMAN_READABLE, LOSS_ATTRIBUTES, GAIN_ATTRIBUTES


if __name__ == "__main__":
    parser = ArgumentParser("Interpretation of Structural Copy Number Variants")

    parser.add_argument("-i", "--input", required=True, help="data to be predicted", type=str)
    parser.add_argument("-c", "--cnv_type", required=True, help="type of cnv, either loss or gain",
                        type=str, choices=["loss", "gain"])
    parser.add_argument("-o", "--output", required=False, help="where the output tsv will be saved",
                        type=str, default="./isv_predictions.tsv")
    parser.add_argument("-p", "--proba", required=False, type=str, choices=["true", "false"],
                        default="true", help="Return probabilities")
    parser.add_argument("-sv", "--shapvalues", type=str, default="false", choices=["true", "false"],
                        help="Calculate SHAP Values")

    args = parser.parse_args()

    if 'tsv' not in args.input:
        sys.exit("Input should be a tsv file")

    if args.input.endswith('.gz'):
        X_raw = pd.read_csv(args.input, sep='\t', compression='gzip')
    elif args.input.endswith('.tsv'):
        X_raw = pd.read_csv(args.input, sep='\t')

    yhat = predict(X_raw, args.cnv_type, proba=(args.proba == "true"))
    res = pd.DataFrame({"ISV": yhat})

    if args.shapvalues == "true":
        attributes = [LOSS_ATTRIBUTES, GAIN_ATTRIBUTES][(args.cnv_type == 'gain') * 1]
        hr_attributes = [HUMAN_READABLE[i] for i in attributes]
        sv = shap_values(X_raw, args.cnv_type)
        res = pd.concat([res, pd.DataFrame(sv.values, columns=hr_attributes)], axis=1)

    res.to_csv(args.output, sep='\t')
