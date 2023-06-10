# trd_root
TRD analysis made in ROOT for beam tests and others

This repository contains files for analysis of flat root file output of [JANA4ML4FPGA](https://github.com/JeffersonLab/JANA4ML4FPGA) software used in various beam tests. 

# Data

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
   
    

Give your token a descriptive name in the "Note" field.
Under "Select scopes", check the box next to "repo" to give the token full control of private repositories.
Click "Generate token".
Copy the token to your clipboard. For security reasons, after you navigate off the page, you will not be able to see the token again.
