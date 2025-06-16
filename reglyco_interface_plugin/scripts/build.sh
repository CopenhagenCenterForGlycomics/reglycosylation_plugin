#!/bin/bash

curl 'https://glyco.me/docs/css/glycosuite.css' > web/glyco_css/glycosuite.css
curl 'https://glyco.me/docs/css/header.css' > web/glyco_css/header.css
curl 'https://glyco.me/docs/css/pageheader_common.css' > web/glyco_css/pageheader_common.css


./node_modules/.bin/webpack --context ./node_modules/sviewer --config ./node_modules/sviewer/webpack.config.js --output web/sviewer.bundle.js
npx lightningcss --bundle --targets "chrome 78" web/main.css  -o  web/main_bundled.css