# diffBase

## Dependencies

- `BLAST+` (capable of generating version 4 BLAST databases)
- `R` version 3.6.x or later

## R package dependencies

### CRAN

`shiny`, `shinythemes`, `ape`, `tidytree`, `magrittr`, `readr`, `gmailr`

### BioConductor (version 10 or later)

`ggtree`, `Biostrings`

## Gettings emails from the community message function

Right now the apps send emails automatically whenever someone submits a message using the community feature. To do this, it uses local GMail keys. For obvious reasons, these are not included in the GitHub repo.

## Configuration

Each app has a `diffBaseConfig.txt` file which can be edited to change some of the behaviours of the apps.

- `URL              = "https://diffbase.uwaterloo.ca"`: server URL
- `UseGmail         = TRUE`: Either `TRUE` or `FALSE`, depending on whether you want email enabled. Useful if you want to run the apps locally and don't want to have to edit the code to prevent it from trying to send emails.
- `ServerEmail      = "diffbaseserver@gmail.com"`: server email
- `GmailCredentials = "~/diff-base/credentials.json"`: GMail credentials. Look up the `gmailr` package for instructions on how to set this up.
- `GmailCache       = "~/diff-base/.secret"`: GMail cache. Look up the `gmailr` R package for instructions on how to set this up.
- `UseBlastp        = TRUE`: Either `TRUE` or `FALSE`, depending on whether you want to allow BLASTP searches.
- `BlastpParameters = "-evalue 1"`: Tune the BLASTP search parameters.

## Running the apps locally

Running the entire app locally isn't possible, but you can run the apps individually with:

```r
shiny::runApp("app-A")
shiny::runApp("app-B")
```

## Running the app with Shiny Server open source edition

Edit your `/etc/shiny-server/shiny-server.conf` (replace `$USER` with your username):

```
run_as shiny;

sanitize_errors false;
preserve_logs true;

server {
  listen 80;
  location / {
      run_as $USER;
      directory_index on;
      app_idle_timeout 0;
      site_dir /home/$USER/diffBase/index/;
      log_dir /home/$USER/diffBase-logs/index/;
    }
  location /app-A {
      run_as $USER;
      app_idle_timeout 0;
      app_dir /home/$USER/diffBase/app-A/;
      log_dir /home/$USER/diffBase-logs/app-A/;
    }
  location /app-B {
      run_as $USER;
      app_idle_timeout 0;
      app_dir /home/$USER/diffBase/app-B/;
      log_dir /home/$USER/diffBase-logs/app-B/;
    }
}
```

