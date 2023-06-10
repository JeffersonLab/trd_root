# trd_root
TRD analysis made in ROOT for beam tests and others

This repository contains files for analysis of flat root file output of [JANA4ML4FPGA](https://github.com/JeffersonLab/JANA4ML4FPGA) software used in various beam tests. 

## Data

- DATA /home/hdtrdops/DATA
- LOG /gluonraid3/data4/rawdata/trd/LOG
- ROOT /home/hdtrdops/ROOT

## Using Git on Gluons

Gluon machines have some nasty firewall settings. You can clone git repos only using https protocol.
Which means that you have to use GitHub repo links starting with `https://...` and use tokens as password
to push your changes. To get your token: 

### Generate a GitHub Personal Access Token (PAT):

First, you need to create a GitHub token that will be used to authenticate your Git operations.

- Use this link to [open GitHub token settings](https://github.com/settings/tokens) or:
    - Log in to your GitHub account.
    - Click on your profile picture in the top right corner and choose "Settings".
    - In the left sidebar, click "Developer settings".
    - Click "Personal access tokens".

- Click "Generate new token"

   ![create token](https://github.com/JeffersonLab/trd_root/blob/main/doc/git_create_token.png?raw=true)

- Give your token a descriptive name in the "Gluon trd_reco" field.
- Under "Select scopes", check the box next to "repo" - "public_repo" to give the token right to push to public repos
- Click "Generate token".
- Copy the token to your clipboard. **SAVE IT! For security reasons, after you navigate off the page, you will not be able to see the token again. **

![create token](https://github.com/JeffersonLab/trd_root/blob/main/doc/git_copy_token.png?raw=true)

*(this token is shown in demonstration purposes and is already delted, so nothing is exposed in this image)

### Use tokens

Now when you do `git push` it will ask:
- username (use your GitHub username)
- password (use your TOKEN here as if it is a password)

Example: 
```bash
[hdtrdops@gluon100 ~/trd_root]$ git push
Username for 'https://github.com': DarTeots
Password for 'https://DarTeots@github.com':
```

### More info

Now Small FAQ:

- [How not to enter username - password every time?](https://git-scm.com/docs/gitcredentials)
- [Is there a safe way?](https://git-scm.com/book/en/v2/Git-Tools-Credential-Storage)
- [What are scopes of auth. for Tokens?](https://docs.github.com/en/apps/oauth-apps/building-oauth-apps/scopes-for-oauth-apps)