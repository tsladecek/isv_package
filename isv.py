from argparse import ArgumentParser
import pandas as pd

from app.scripts.predict import Predict


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

    p = Predict(cnv_type=args.cnv_type,
                data_path=args.input)

    yhat = p.predict(proba=args.proba)
    sv = p.shap_values().values

    print(sv)
    # pd.DataFrame({"ISV": yhat}).to_csv(args.output, sep='\t')
