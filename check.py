from pathlib import Path

INPUT_DIR = Path(
    r"C:\Users\bhatt\Desktop\Omics\diversification_barley\Assemblytics"
)

files = sorted(
    INPUT_DIR.rglob("*Assemblytics_structural_variants.bed")
)

print("=" * 80)
print("ASSEMBLYTICS FILE INSPECTION")
print("=" * 80)
print(f"Folder: {INPUT_DIR}")
print(f"Files found: {len(files)}")

for i, file in enumerate(files[:5], 1):

    print("\n" + "-" * 80)
    print(f"FILE {i}: {file.name}")
    print("-" * 80)

    try:
        with open(
            file,
            "r",
            encoding="utf-8-sig",
            errors="replace"
        ) as f:

            for line_number in range(1, 6):

                line = f.readline()

                if not line:
                    break

                print(
                    f"LINE {line_number}: "
                    f"{repr(line[:1000])}"
                )

    except Exception as e:
        print("ERROR:", e)

print("\nInspection finished.")