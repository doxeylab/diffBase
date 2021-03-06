# diffBase

## Dependencies

- `BLAST+` (capable of generating version 4 BLAST databases)
- `R` version 3.6.x or later

## R package dependencies

### CRAN

`shiny`, `shinythemes`, `ape`, `tidytree`, `magrittr`, `readr`, `gmailr`

### BioConductor (version 10 or later)

`ggtree`, `Biostrings`

## Updating the database

The update script `update.R` can be run to refresh the data files in the apps. The script can handle changes in the sequences themselves as well as the number of subtype sequences. It **cannot** handle changes in the number of group families. For these cases, you must manually update the source code of the apps and the update script.

The script checks for changes in the sequences by comparing them to a store of MD5 checksums, and only executes if it finds any changes. A couple of environmental variables can be set:

- `FORCE_UPDATE=T`: Force the update script to run, even if the MD5 checksums haven't changed. Useful if you've changed the behaviour of the update script or if you've updated the metadata.
- `NO_UPDATE_TREE=T`: Don't update the trees. Useful if you manually generated the trees and don't want to build them using the mechanism in the update script.

Note that the update script requires the `phangorn` R package from CRAN to build the trees.

Additionally, note that the subtypes are named based on their order within the family files. So if you wish to add subtypes and not disturb the current naming scheme, add them to the bottom of the fasta files.

## Gettings emails from the community message function

Right now the apps send emails automatically whenever someone submits a message using the community feature. To do this, it uses local GMail keys. For obvious reasons, these are not included in the GitHub repo.

## Configuration

Each app has a `diffBaseConfig.txt` file which can be edited to change some of the behaviours of the apps. Remember to change the file paths to match your system.

- `URL              = "https://diffbase.uwaterloo.ca"`: server URL
- `UseGmail         = TRUE`: Either `TRUE` or `FALSE`, depending on whether you want email enabled. Useful if you want to run the apps locally and don't want to have to edit the code to prevent it from trying to send emails.
- `ServerEmail      = "diffbaseserver@gmail.com"`: server email
- `GmailCredentials = "/home/$USER/diffBase/credentials.json"`: GMail credentials. Look up the [`gmailr`](https://github.com/r-lib/gmailr) package for instructions on how to set this up.
- `GmailCache       = "/home/$USER/diffBase/.secret"`: GMail cache. Look up the [`gmailr`](https://github.com/r-lib/gmailr) R package for instructions on how to set this up.
- `UseBlastp        = TRUE`: Either `TRUE` or `FALSE`, depending on whether you want to allow BLASTP searches.
- `BlastpParameters = "-evalue 1"`: Tune the BLASTP search parameters.

## Running the apps locally

Running the entire app locally isn't possible, but you can run the apps individually with:

```r
shiny::runApp("app-A")
shiny::runApp("app-B")
```

## Running the app with Shiny Server open source edition

Note: these instructions were written with Ubuntu 16.04.6 LTS in mind. They also assume you've created the log folders.

### No SSL

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

If you ever make changes to the app, restart the server:

```
sudo systemctl restart shiny-server
```

### With SSL

If you've installed Shiny Server already, stop it:

```
sudo systemctl stop shiny-server
```

Install `nginx` and stop it:

```
sudo apt-get install nginx
sudo systemctl stop nginx
```

Now edit their configs:

```
# /etc/nginx/sites-available/default
server {
  listen 80;
  return 301 https://$host$request_uri;
}
server {
  listen 443 ssl;
  server_name diffbase.uwaterloo.ca;
  ssl_certificate /etc/letsencrypt/live/diffbase.uwaterloo.ca/fullchain.pem;
  ssl_certificate_key /etc/letsencrypt/live/diffbase.uwaterloo.ca/privkey.pem;
  location / {
    proxy_pass http://127.0.0.1:4949;
  }
}
```

Replace `$USER` with the server username:

```
# /etc/shiny-server/shiny-server.conf
run_as shiny;
sanitize_errors false;
preserve_logs true;
server {
  listen 4949 127.0.0.1;
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

Install `certbot` (based on these [instructions](https://geekflare.com/setup-nginx-with-lets-encrypt-cert/)):

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:certbot/certbot
sudo apt-get update
sudo apt-get install python-certbot-nginx
```

Run `certbot`:

```
sudo certbot --nginx
```

It will ask for a domain name. For example, `diffbase.uwaterloo.ca`. Once it's done it will tell you when the certificate will expire. When the time comes, renew it with:

```
sudo certbot renew
```

The `certbot` program will try and automatically modify your `nginx` config file. Make sure it did it correctly.

Start `nginx` and `shiny-server`:

```
sudo systemctl start nginx
sudo systemctl start shiny-server
```
