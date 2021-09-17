import markdown

with open("README.md", "r", encoding="utf-8") as input_file:
    text = input_file.read()
html = markdown.markdown(text, extensions=['tables','nl2br','extra'])
print (html)
with open("some_file.html", "w", encoding="utf-8", errors="xmlcharrefreplace") as output_file:
    output_file.write(html)