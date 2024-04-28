# 博士毕业论文 LaTeX 源码

## 文件管理方式

- 主要宏包为 `marco.tex`，它包含所有 `documents` 即正文环境以外的所有定义；
- `thesis.bib` 是所有文献的打包；
- 一级目录结构：
  - `chap-*` 是管理正文章节的目录；
  - `dummy` 对每一章节生成一份 pdf 文件；
  - `intro` 是管理摘要、后记的目录；
  - `notations` 是管理附录的软件、方法列表，以及术语对照表的目录；
  - `draft` 是统合所有分散的章节，统稿为完整论文的目录。
- 二级目录结构：
  - `dummy-*.tex` (正文章节在文件夹 `dummy` 中、其他补充部分在各子目录中) 是每个子章节 (或附录) 的外壳文件。每个章节可以独立编译出 pdf 文件，这些 `dummy-*.tex` 是实际被编译的源码；但正文的具体内容不在这些外壳文件中。
  - `chap-*.tex`、`appendix-*.tex` 或其他 .tex 文件是正文文件。所有重要的内容都在这些文件中；但它们不应该包含在 `\begin{document}` 环境中，而是用其他文件作为外壳。
  - `assets` 是图像文件、以及生成图像的源代码和原始数据文件。

## 工作流

本文档编译工作流是 VSCode + Latex Workshop。合理配置的话，应该在保存正文文件时，可以自动生成**对应章节**的 pdf 文档 (不是全文的文档)。

全文文档需要编译 `draft/draft.tex` 文件；这个需要手动触发。编译该文档在高速 NVME SSD 上通常需要半分钟到一分钟。

VSCode 工作需要注意的事项包括
- 更改默认编译器为 xelatexmk；
- 编译间间隔不要太短，不然刚刚打完一两个词，笔记本的风扇就会狂转；
- 保存文件的时候就编译了，不用刻意到命令行界面；
- 使用内置 pdf 浏览器时，Ctrl + Left Mouse 可以 pdf -> tex；Ctrl + Alt + J 可以 tex -> pdf；
- 在正文 .tex 文件第一行写上对应的外壳文件 (比如第一章是 `% !TEX root=../dummy/dummy-01.tex`)；这会告诉 Latex Workshop 应该要编译哪一个文件。

建议使用 TexLive 2024 编译，或将 `fduthesis` 升级到 v0.9a。否则 `siunitx` 和 `dcolumn` 包的存在可能会干扰编译。

对于文献管理，
- 依照现在的 `gbt7714` 库，在条目中有 doi 项目时，请尽量不要引入 url 与 urldata 两栏；
- 参考文献需要校对；请特别注意校对标题是否显示正常，以及含有非 ASCII 拉丁字母的姓名是否正常地大写了；
- 可能 biblatex 比 bibtex 要更方便，但目前直接使用 biblatex 会报错。

若涉及到跨章节交叉引用，可以考虑在 `marco.tex` 与 `draft.tex` 中，定义与重定义的 `\alertref` 命令。该命令在编译独立章节文件中，只是标红的字符串；在编译完整文档时，该命令是交叉引用。

## 附录

VSCode 用户设定

```json
{
    "workbench.iconTheme": "catppuccin-latte",
    "workbench.colorTheme": "Catppuccin Latte",
    "editor.fontFamily": "Sarasa Fixed SC, Consolas, 'Courier New', monospace",
    "latex-workshop.latex.recipes": [

        {
            "name": "latexmk (xelatex)",
            "tools": [
                "xelatexmk"
            ]
        },
        {
            "name": "latexmk",
            "tools": [
                "latexmk"
            ]
        },
        {
            "name": "latexmk (latexmkrc)",
            "tools": [
                "latexmk_rconly"
            ]
        },
        {
            "name": "latexmk (lualatex)",
            "tools": [
                "lualatexmk"
            ]
        },
        {
            "name": "pdflatex -> bibtex -> pdflatex * 2",
            "tools": [
                "pdflatex",
                "bibtex",
                "pdflatex",
                "pdflatex"
            ]
        },
        {
            "name": "Compile Rnw files",
            "tools": [
                "rnw2tex",
                "latexmk"
            ]
        },
        {
            "name": "Compile Jnw files",
            "tools": [
                "jnw2tex",
                "latexmk"
            ]
        },
        {
            "name": "Compile Pnw files",
            "tools": [
                "pnw2tex",
                "latexmk"
            ]
        },
        {
            "name": "tectonic",
            "tools": [
                "tectonic"
            ]
        }
    ],
    "editor.wordWrap": "on",
    "editor.unicodeHighlight.allowedLocales": {
        "zh-hant": true,
        "zh-hans": true
    },
    "latex-workshop.view.pdf.external.viewer.command": "C:\\Program Files\\SumatraPDF\\SumatraPDF.exe",
    "latex-workshop.view.pdf.external.synctex.command": "C:\\Program Files\\SumatraPDF\\SumatraPDF.exe",
    "editor.fontSize": 15,
    "editor.cursorBlinking": "smooth",
    "editor.cursorSmoothCaretAnimation": "explicit",
    "workbench.list.smoothScrolling": true,
    "files.autoSave": "afterDelay",
    "files.autoSaveDelay": 10000,
    "latex-workshop.latex.autoBuild.interval": 5000,
    "latex-workshop.latex.autoBuild.run": "onSave",
    "security.workspace.trust.untrustedFiles": "open",
    "latex-workshop.bibtex-format.tab": "4",
    "security.allowedUNCHosts": [
        "wsl.localhost"
    ],
    "editor.minimap.renderCharacters": false,
}
```
