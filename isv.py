from argparse import ArgumentParser
import pandas as pd

from app.scripts.predict import predict


if __name__ == "__main__":
    parser = ArgumentParser("Interpretation of Structural Copy Number Variants")

    parser.add_argument("-i", "--input", required=True, help="data to be predicted", type=str)
    parser.add_argument("-c", "--cnv_type", required=True, help="type of cnv, either loss or gain",
                        type=str, choices=["loss", "gain"])
    parser.add_argument("-o", "--output", required=False, help="where the output tsv will be saved",
                        type=str, default="./isv_predictions.tsv")

    args = parser.parse_args()

    y = predict(cnv_type=args.cnv_type,
                data_path=args.input)

    pd.DataFrame({"ISV": y}).to_csv(args.output, sep='\t')
