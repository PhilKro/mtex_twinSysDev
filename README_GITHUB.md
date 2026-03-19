# MTEX Custom Development Environment

This repository contains a custom MTEX development workspace. The environment is specifically architected to isolate code from heavy data, preventing Windows file locks and keeping cloud syncs stable.

## 1. System Architecture

This setup relies on physical separation and a custom Git routing strategy.

* **The Local Workspace:** All code execution and version control happens in the local `C:\Users\phkr\MTEX_workEnv` folder. This directory is strictly local and completely outside of OneDrive.

* **The Data Junction:** A Windows Directory Junction named `onedriveData` acts as a portal inside the workspace. It looks like a local folder to MATLAB but physically reads and writes directly to `C:\Users\phkr\OneDrive - empa.ch\MTEX_Data`.

* **The `.gitignore`:** A heavily customized ignore file ensures that Git ignores the `onedriveData` portal, all heavy data files (`.ang`, `.ctf`, `.hdf5`), and your custom MATLAB `.mat` files.

### Remote Connections

The local repository communicates with three distinct remote locations:

1. `upstream`: The official MTEX toolbox repository (Used strictly for pulling updates).

2. `origin`: The public GitHub fork (`PhilKro/mtex_twinSysDev`).

3. `onedrive`: A private, "bare" Git repository hidden inside OneDrive, functioning as a daily cloud backup.

## 2. Environment Setup & Recreation

If this environment ever needs to be recreated on a new machine, follow these steps.

**Step 1: Clone and Branch**

    git clone https://github.com/PhilKro/mtex_twinSysDev.git "C:\Users\phkr\MTEX_workEnv"
    cd "C:\Users\phkr\MTEX_workEnv"
    git checkout -b mtex_workEnv

**Step 2: Link Remotes**

    git remote add upstream https://github.com/mtex-toolbox/mtex.git
    git remote add onedrive "C:\Users\phkr\OneDrive - empa.ch\mtex_private_backup.git"
    git push -u onedrive mtex_workEnv

**Step 3: Create the Data Junction (Run in PowerShell)**

    New-Item -ItemType Junction -Path "onedriveData" -Target "C:\Users\phkr\OneDrive - empa.ch\MTEX_Data"

## 3. Daily Workflow

For day-to-day coding, work entirely inside the local `MTEX_workEnv` folder.

### Accessing Data in MATLAB

When writing scripts, use the relative path through the junction portal. MATLAB will automatically reach into OneDrive:

    % Example data load using the junction portal
    ebsd = EBSD.load('onedriveData/my_scan.ang', 'interface', 'ang');

### Saving Your Work

The push command is hardwired to send your code to the private OneDrive backup. When you finish a session, run these commands in the terminal (or use the VS Code Sync button):

    git add .
    git commit -m "Describe your changes"
    git push

## 4. Updating MTEX (The Rebase Workflow)

When the official MTEX developers release new features, use a rebase to update the underlying code while preserving your custom twin system edits on top.

**Step 1: Fetch the latest official code**

    git fetch upstream

**Step 2: Rebase your branch**
Target the active development branch for the latest features, or a specific version tag for stability:

    # To get the absolute latest features:
    git rebase upstream/develop
    
    # OR to target a specific stable release:
    git rebase tags/mtex-6.1.1

**Step 3: Force push to your backup**
Because a rebase rewrites your local Git history, you must force the update onto your OneDrive backup:

    git push -f onedrive mtex_workEnv

## 5. Diagnostics & Troubleshooting

Use these commands if you need to verify your Git routing.

**Check where "Push" is pointing:**

    git branch -vv

*Expected Output: `* mtex_workEnv [onedrive/mtex_workEnv]`*

**Check all remote URLs:**

    git remote -v

*Expected Output: A list showing `onedrive`, `origin`, and `upstream` with their correct paths.*

### Common Fixes

* **Missing Data:** If MATLAB cannot find your files, the `onedriveData` junction may have been deleted. Recreate it in PowerShell:
  `New-Item -ItemType Junction -Path "onedriveData" -Target "C:\Users\phkr\OneDrive - empa.ch\MTEX_Data"`

* **Merge Conflicts:** If Git pauses during a rebase, MTEX updated a file that you also edited. Open the conflicted file, choose which code to keep, save the file, run `git add <filename>`, and then run `git rebase --continue`.

* **Stuck Rebase:** If Git reports a rebase is already in progress, clear it by running `git rebase --abort`.