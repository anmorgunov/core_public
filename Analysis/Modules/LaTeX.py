from pathlib import Path
from typing import List, Union, cast


class Table:

    @staticmethod
    def export_table(
        caption: str,
        label: str,
        positioning: str,
        headers: Union[str, List[str]],
        body: List[List[str]],
        save_path: str,
    ) -> None:
        file_content = "\\begin{table}\n"
        file_content += "  \caption{\\textbf{%s}}\n" % caption
        file_content += "  \label{tbl:%s}\n" % label
        file_content += "  \\begin{tabular}{%s}\n" % positioning
        file_content += "    \\toprule\n"

        def printrow(row: List[str]) -> List[str]:
            return ["    "] + [
                f"{rowElt} & " if i + 1 != len(row) else f"{rowElt} \\\ \n"
                for i, rowElt in enumerate(row)
            ]

        if isinstance(headers[0], list):
            for header in headers:
                file_content += "".join(printrow(header))
            assert len(positioning.split()) == len(
                header
            ), "You must specify the alignment of each column"
        else:
            if "|" in positioning:
                splitted = positioning.split(" | ")
            else:
                splitted = positioning.split()
            assert len(splitted) == len(
                headers
            ), "You must specify the alignment of each column"
            file_content += "".join(printrow(cast(List[str], headers)))
        file_content += "    \midrule\n"
        for row in body:
            file_content += "".join(printrow(row))
        file_content += "    \\bottomrule\n"
        file_content += "  \end{tabular}\n"
        file_content += "\end{table}\n"
        with open(save_path, "w") as f:
            f.write(file_content)


if __name__ == "__main__":
    testHead = ["Molecule", "Atom", "UHF", "UMP2"]
    testBody = [["co", "o", "12.3", "12.1"], ["co", "c", "8.5", "8.4"]]
    testCaption = "This is an example of a table"
    testLabel = "sampletable"
    testPositioning = "c c c c"
    testPath = Path(__file__).resolve().parent / "testTable.tex"
    Table.export_table(
        testCaption, testLabel, testPositioning, testHead, testBody, str(testPath)
    )
