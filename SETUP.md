# Project Template for RStudio Server Docker container with renv and tidyverse

## Steps for Geting Started

- Run the `nf-core` pipeline on your raw data. The `nf-core` pipeline is a community-curated collection of analysis pipelines built using Nextflow, designed for reproducibility and best practices in bioinformatics. For more details and installation instructions, visit [nf-core documentation](https://nf-co.re/).
- Clone the repo to your project directory
`git clone https://<PERSONALACCESSTOKENS>@github.com/CECADBioinformaticsCoreFacility/ProjectTemplate.git`

> **Note:** Replace `<PERSONALACCESSTOKENS>` with your actual GitHub personal access token. Ensure that your personal access token is kept private and not shared or exposed in public repositories. You can generate a personal access token by following the instructions in the [GitHub documentation](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token).

- Customise the `compose.yml` file to set the environment variables and service configurations as per your project requirements.
    - `environment:`
    - `volumes:`
> **Note:** You can check userid and groupid by `id -u` and `id -g` in the terminal

- The working directory gets mounted in the guest system as `/home/rstudio/project` and the `${HOME}/.cache` directory gets mounted as `/home/rstudio/.cache` by default 
- Use `docker compose up -d` to start the container
- use `docker compose down` to stop the running container
- Use `docker compose up -d --build` rebuild the container and if you modefy the Docker file inside `docker/` directory
- add your user to the docker group instead of running Docker with sudo

```
# check /etc/group to see if the docker group exists
cat /etc/group | grep docker

# create a docker group if it does not exist
sudo groupadd docker

# add yourself to the docker group
sudo usermod -aG docker $(whoami)
```


