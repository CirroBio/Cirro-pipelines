from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def main():

    # Instantiate the Cirro dataset object
    ds = PreprocessDataset.from_running()

    # Set up a CSV with the regions provided as inputs
    regions = find_regions(ds)

    # Save to a file
    regions.to_csv("regions.csv", index=None)

    # Log the contents of the CSV
    ds.logger.info("Regions:")
    for line in regions.to_csv(index=None).split("\n"):
        ds.logger.info(line)


def find_regions(ds: PreprocessDataset) -> pd.DataFrame:
    return pd.DataFrame([
        dict(
            id=input["name"],
            uri=input["s3"] + "/data/region.json"
        )
        for input in ds.metadata["inputs"]
    ])


if __name__ == "__main__":
    main()
