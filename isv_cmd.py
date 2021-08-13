import argparse
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import sys

from isv.annotate import annotate
from isv.predict import predict
from isv.shap_vals import shap_values
from isv.scripts.constants import HUMAN_READABLE, LOSS_ATTRIBUTES, GAIN_ATTRIBUTES


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

    loss_indices = np.where(bed.iloc[:, 3] == "DEL")[0]
    gain_indices = np.where(bed.iloc[:, 3] == "DUP")[0]

    assert len(loss_indices) + len(gain_indices) == len(bed),\
        "All CNVs should have cnv type of either 'DUP' or 'DEL'"

    annotated = annotate(bed)

    final_predictions = np.empty(len(bed), dtype=[np.int8, np.float64][args.proba * 1])

    print("Predicting")
    if loss_indices.shape[0] > 0:
        final_predictions[loss_indices] = predict(annotated.iloc[loss_indices], "loss", proba=args.proba)
    if gain_indices.shape[0] > 0:
        final_predictions[gain_indices] = predict(annotated.iloc[gain_indices], "gain", proba=args.proba)

    bed["ISV"] = final_predictions

    if args.shapvalues:
        print("Calculating SHAP values")
        res = []

        if loss_indices.shape[0] > 0:
            attributes = LOSS_ATTRIBUTES
            hr_attributes = ['SHAP_' + HUMAN_READABLE[i].replace(' ', '_') for i in attributes]
            sv = shap_values(annotated.iloc[loss_indices], "loss")
            res.append(pd.DataFrame(sv.values, columns=hr_attributes))

        if gain_indices.shape[0] > 0:
            attributes = GAIN_ATTRIBUTES
            hr_attributes = ['SHAP_' + HUMAN_READABLE[i].replace(' ', '_') for i in attributes]
            sv = shap_values(annotated.iloc[gain_indices], "gain")
            res.append(pd.DataFrame(sv.values, columns=hr_attributes))

        res = pd.concat(res)

        res.index = np.concatenate([loss_indices, gain_indices])

        res = res.sort_index()
        # Append to results
        bed = pd.concat([bed, res], axis=1)

    bed.to_csv(args.output, sep='\t', index=False)
    print(f"Results saved to {args.output}")
