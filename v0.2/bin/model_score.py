import pandas as pd
import joblib
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Print to stdout for Nextflow logs
    print(">>> Running model_score.py for", args.output, flush=True)

    debug_messages = []

    try:
        debug_messages.append(f"Loading features from: {args.features}")
        X = pd.read_csv(args.features, sep="\t", index_col=0)
        debug_messages.append(f"Features shape: {X.shape}")

        debug_messages.append(f"Loading model from: {args.model}")
        model_dict = joblib.load(args.model)
        scaler = model_dict["scaler"]
        model = model_dict["model"]
        feature_names = model_dict["feature_names"]
        debug_messages.append(f"Model features: {feature_names}")

        missing = set(feature_names) - set(X.columns)
        if missing:
            raise ValueError(f"Missing features in input data: {missing}")

        X_ordered = X[feature_names]
        debug_messages.append(f"Ordered features shape: {X_ordered.shape}")

        X_scaled = scaler.transform(X_ordered)
        debug_messages.append(f"Scaled features shape: {X_scaled.shape}")

        y_pred = model.predict(X_scaled)
        debug_messages.append(f"Predicted values shape: {y_pred.shape}")

        output_df = pd.DataFrame({
            "sample_id": X.index,
            "predicted_expression": y_pred
        }).set_index("sample_id")

        output_df.to_csv(args.output, sep="\t")
        debug_messages.append(f"Saved predicted scores to: {args.output}")

    except Exception as e:
        debug_messages.append(f"Error: {e}")

    debug_file = args.output + ".debug.txt"
    with open(debug_file, "w") as f:
        for msg in debug_messages:
            f.write(msg + "\n")

if __name__ == "__main__":
    main()
