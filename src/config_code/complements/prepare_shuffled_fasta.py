import random
from pathlib import Path

from Bio import SeqIO


PROJECT_PATH = Path("C:/Users/chait/Downloads/BINF/BINF_2111")

VIRAL_PATH = (
    PROJECT_PATH
    / "src/config_code/hpv_step2/targets"
    / "HPV16_E6starI_E7_spliced_region.rna.fasta"
)

OUTPUT_DIR = PROJECT_PATH / "results/step5/shuffled_controls"

NUM_SHUFFLES = 100
RANDOM_SEED = 42
BATCH_SIZE = 50


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    record = SeqIO.read(VIRAL_PATH, "fasta")
    viral_sequence = str(record.seq).upper().replace("T", "U")

    rng = random.Random(RANDOM_SEED)

    for batch_start in range(0, NUM_SHUFFLES, BATCH_SIZE):
        batch_number = batch_start // BATCH_SIZE + 1
        output_path = OUTPUT_DIR / f"viral_shuffles_batch_{batch_number}.fasta"

        with output_path.open("w", encoding="utf-8") as output:
            batch_end = min(batch_start + BATCH_SIZE, NUM_SHUFFLES)

            for index in range(batch_start, batch_end):
                bases = list(viral_sequence)
                rng.shuffle(bases)
                shuffled_sequence = "".join(bases)

                output.write(
                    f">{record.id}_shuffle_{index + 1:03d}\n"
                    f"{shuffled_sequence}\n"
                )

        print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()