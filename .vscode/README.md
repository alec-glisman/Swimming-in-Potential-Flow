# Directory: .vscode

Files relevant for developing the project in VSCode.  

## Files

`c_cpp_properties.json`: Provides VSCode information on where to find header files, compilers, and standards used in compiling the project.
Current have a setup for both MacOS and Linux.  

`extensions.json`: List of all extensions used in development.
This incldes both necessary extensions, such as C++ intellisense, and quality of life extensions, such as color themes.

`launch.json`: Configurations for debugging files.

`settings.json`: VSCode settings for this project  

`tasks.json`: Tasks to call in VSCode.
These build (and run) the project in a range of configurations.

## Useful Keyboard Shortcuts

### Navigation

#### Command Pallete

- `^⇧P` Open command pallette (`>` command)
- `^P` Search for files in current workspace
- `^T` Search for symbols in current workspace (`@` command)
- `^⇧O` Search for symbols in current file (`#` command)

#### Definition/References

- `F12` Go to definition
- `⌥F12` Peek definition
- `⇧F12` Peek references (`⌘k` and `F12` on MacOS)
- `^Tab` Switch between open tabs in current editor
- `^<number>` Swtich to editor group number (`0` is sidebar)

### Auto-completion

- `^⎵` Autocomplete suggestions
- `^⇧⎵` Function assist (shows function signature)

### Editing

- `F2` Change name of symbol
- `F8` See intellisense errors
- `^.` Show recommended actions for currently highlighted code
- `⌥⇧ ←/→` Select blocks of code (`⌘^ ⇧ ←/→`)
- `^⇧\` Navigate to matching bracket
- `^⌥↓`/`^⌥↑` Add multiple cursors
- `^\` Split editor

### Terminal

- `^`\` Open/close terminal pane
- `^⇧`\` Open new terminal instance
