# CHANGELOG



## v0.4.3 (2025-08-20)

### Chore

* chore: update poetry lock ([`3fbd379`](https://github.com/lanl/ThunderBoltz/commit/3fbd379014ab962745625174ae86e9e04df6802c))

### Fix

* fix: Limit mpl to avoid further deprecations ([`91213ae`](https://github.com/lanl/ThunderBoltz/commit/91213ae587ed4ca4c040e3c24f97840b98db4668))


## v0.4.2 (2025-08-20)

### Documentation

* docs: update README with article

and fix the example with out of date directory. ([`21cf69c`](https://github.com/lanl/ThunderBoltz/commit/21cf69cf1dc331476318d97d6cd6967e1470bdb9))

### Fix

* fix: Remove redundant import for python 3.10 ([`ed71d11`](https://github.com/lanl/ThunderBoltz/commit/ed71d118e1fc2dd4b46108d597c6f70376ad3502))

### Unknown

* doc: update diffusion docs in api manual ([`6fe8a12`](https://github.com/lanl/ThunderBoltz/commit/6fe8a124aa1375d329804cbdc00e5150b3990891))


## v0.4.1 (2024-07-19)

### Fix

* fix: use c++ rng instead of c built-in ([`668ee6f`](https://github.com/lanl/ThunderBoltz/commit/668ee6f452a2885e7a523d1524e64347ffa1974b))

* fix: Nmin criteria is now taken as max

Instead of using N_A N_B product criteria to determine NP,
it will only do so if more particles are required compared to the NP
set by the user. ([`8056016`](https://github.com/lanl/ThunderBoltz/commit/8056016ee868cdbfaa38c4fe094e0d7e6530f13f))


## v0.4.0 (2024-07-10)

### Chore

* chore: forward difference takes default parameter ([`564cd2b`](https://github.com/lanl/ThunderBoltz/commit/564cd2b28ae8b5eeed0fc16133de9793d09d3a74))

### Documentation

* docs: remove reference to linear model

From old curve fit code (didn&#39;t work because of small time step
values). ([`2420449`](https://github.com/lanl/ThunderBoltz/commit/242044990634961f150c5cc84a66cd4458d50e61))

### Feature

* feat: Add user defined timeseries plotting. ([`862ec81`](https://github.com/lanl/ThunderBoltz/commit/862ec81e8c94b767225a83427c8fc0c15c3016c3))

* feat: added plot_tree ([`bd9534d`](https://github.com/lanl/ThunderBoltz/commit/bd9534d239752b23523c773f8b136c1c12e33ab2))

### Test

* test: update test for fits -&gt; fit rename ([`7bc391c`](https://github.com/lanl/ThunderBoltz/commit/7bc391c9c08382a2f2c41bf1669f1697093ba2cb))


## v0.3.0 (2024-06-22)

### Documentation

* docs: update readme + spelling ([`656eeb4`](https://github.com/lanl/ThunderBoltz/commit/656eeb464e75ba81df805974fff91be152a3a8a0))

* docs: spelling ([`94fbb13`](https://github.com/lanl/ThunderBoltz/commit/94fbb130037977ad4d917406ecf042b8294205ed))

* docs: update onsager analytic rate expression ([`4fa850f`](https://github.com/lanl/ThunderBoltz/commit/4fa850f8b3c2b1d69847ff428cb7f7c1ef24c07f))

* docs: add diffusion documentation ([`9ff3b33`](https://github.com/lanl/ThunderBoltz/commit/9ff3b33fcb2626215ac7adde483671e61fd5e4d4))

### Feature

* feat: allow cross section plotter to accept labels ([`1e48e80`](https://github.com/lanl/ThunderBoltz/commit/1e48e803d64146d14c2f4af96130f3ae5eb71f6e))

* feat: add optional fit for all time derivatives ([`644afc2`](https://github.com/lanl/ThunderBoltz/commit/644afc2a5ba63016266be280ddf6fca950fd8d59))

* feat: add Reid ramp setup through python API ([`76b4714`](https://github.com/lanl/ThunderBoltz/commit/76b4714c786fb63c74c9fc83bf730feb1bc2a07b))

* feat: add diffusion correlation calcs to c++ src ([`0049912`](https://github.com/lanl/ThunderBoltz/commit/00499122be5e670f0fdc8184ee007acaac1070ef))

### Fix

* fix: Remove SST experiments from transport figure ([`a7573b6`](https://github.com/lanl/ThunderBoltz/commit/a7573b641df08c86dc4a8e2c7db56bd7ae3043c0))

* fix: fillna preventing ionization calculation ([`4e968ff`](https://github.com/lanl/ThunderBoltz/commit/4e968ff4de0f3031070126bd09bf89f206d01dd9))

* fix: change onsager cross sections and fig ([`924fce9`](https://github.com/lanl/ThunderBoltz/commit/924fce92bb1bc4e6b6ec33598d644f3ac0241832))

* fix: change CR default to 0 (match w/ cpp default) ([`efec58f`](https://github.com/lanl/ThunderBoltz/commit/efec58f5777f042da1fa019607901e1c9106f92a))

### Test

* test: add ExB direction check ([`be8843e`](https://github.com/lanl/ThunderBoltz/commit/be8843ed3121e1b928769c2b146b5866176a4d04))


## v0.2.0 (2024-01-07)

### Documentation

* docs: add open option for doc build file ([`5b35206`](https://github.com/lanl/ThunderBoltz/commit/5b3520626740bd47c55f1f3e4caf8552e6e720c9))

* docs: add more references to cumulative counts ([`a7df529`](https://github.com/lanl/ThunderBoltz/commit/a7df529b87f899d339f3739439c7aa2b0a7584a1))

### Feature

* feat: add method get_counts() and docs ([`a23f110`](https://github.com/lanl/ThunderBoltz/commit/a23f11022d3e230deaf1f6db0e845561ed6e2a58))


## v0.1.2 (2023-12-03)

### Documentation

* docs: add figures to quick start guide ([`e9a8cdc`](https://github.com/lanl/ThunderBoltz/commit/e9a8cdcf7a27342dcbc5ac36eb99e7980391e888))

* docs: fix readme rst in md, add doc badge ([`1aded9b`](https://github.com/lanl/ThunderBoltz/commit/1aded9b12395915f45b98a9dbd3c1a510b92f36f))

* docs: update docs with new repo location

Also, add status badges to readme. ([`89f76c2`](https://github.com/lanl/ThunderBoltz/commit/89f76c20a9c79281ab7647ef39ad6b35444c5464))

* docs: update html short index with pip inst. ([`290decd`](https://github.com/lanl/ThunderBoltz/commit/290decd0ac0d9b6b020362ce5c851925a4b1c85e))

* docs: update html index with pip install ([`72b1faa`](https://github.com/lanl/ThunderBoltz/commit/72b1faa13acf636390bb06a1f2480777232396a4))

### Fix

* fix: add save option to rate figure ([`273cee4`](https://github.com/lanl/ThunderBoltz/commit/273cee48821c71ada0bdc22d2e925095be21819c))


## v0.1.1 (2023-11-07)

### Build

* build: fix typo in pyproject toml config

cc ([`901d626`](https://github.com/lanl/ThunderBoltz/commit/901d6263e87a92cf089bdc989ae8a73fa0cc6e83))

* build: debug release with verbose ([`f1a857f`](https://github.com/lanl/ThunderBoltz/commit/f1a857f79b25faa9188cc428bbcc4ffc1369d33d))

* build: adjust python scipy requirements ([`61aba35`](https://github.com/lanl/ThunderBoltz/commit/61aba35248296f0ea9fdc85f81011ba707b2eb0c))

* build: use compatible python version in RTD build ([`4e1df9e`](https://github.com/lanl/ThunderBoltz/commit/4e1df9ef507e2c0594c1db3fea5c14ea8242ffd2))

* build: remove release spec from sphinx

Add autoapi dir spec. ([`ca4bc1d`](https://github.com/lanl/ThunderBoltz/commit/ca4bc1d177e87463fcec373b348eeac9f9e276c8))

* build: require thunderbolz in rtd reqs ([`fea855a`](https://github.com/lanl/ThunderBoltz/commit/fea855a39962db92e6fdb3a1c1de685882fcaf8b))

### Documentation

* docs: format readme ([`7d6b6b6`](https://github.com/lanl/ThunderBoltz/commit/7d6b6b6ad4d22b1b6dbdd7e823ebd3a5c75c3670))

* docs: require napoleon build vsn ([`3a0f20f`](https://github.com/lanl/ThunderBoltz/commit/3a0f20f94f44f03568d9330e8803d46d016b7f65))

* docs: require specific doc build vsn in rtd req ([`5d4f953`](https://github.com/lanl/ThunderBoltz/commit/5d4f953363d49d6f258d7cb92d8626b112279fd8))

* docs: correct cpp source file loc in readme ([`409ae0a`](https://github.com/lanl/ThunderBoltz/commit/409ae0a47b9de83c6d2566f128b8223782d90cea))

### Fix

* fix: Update readme with pip install option ([`9fd94bc`](https://github.com/lanl/ThunderBoltz/commit/9fd94bc16e5d824ee0ecb34bb3cd0d68983ba701))


## v0.1.0 (2023-11-07)

### Feature

* feat(license): Correct license with version ([`6d0646f`](https://github.com/lanl/ThunderBoltz/commit/6d0646f29de29f4a9dbfec4acb12eec32e97f592))


## v0.0.1 (2023-11-07)

### Build

* build: add thunderboltz to docs reqs ([`4d6c405`](https://github.com/lanl/ThunderBoltz/commit/4d6c4057c89efcd5cc926bcd1e7ec943997140ad))

* build: require sphinx in text ([`ae00de4`](https://github.com/lanl/ThunderBoltz/commit/ae00de471135e2f3465f0d9a8096b918c31bfd4e))

* build: fix release condition ([`6333e8d`](https://github.com/lanl/ThunderBoltz/commit/6333e8d5ab1722fdc0bb53ee3b56ca221ea0a8ea))

* build: require release for CD step ([`03b4935`](https://github.com/lanl/ThunderBoltz/commit/03b4935cba040fb41f8645bc2c8171ad189ea225))

### Documentation

* docs: try version specific sphinx requirements ([`9f18cb7`](https://github.com/lanl/ThunderBoltz/commit/9f18cb7eb4bc097aadcdaaf035146314ab9ceba9))

* docs: add requirements link in yml ([`d6befea`](https://github.com/lanl/ThunderBoltz/commit/d6befeacd16da64d3f8de50161cd194d7bc3fb2f))

* docs: add pip install in readme ([`ae25d66`](https://github.com/lanl/ThunderBoltz/commit/ae25d66bc350a4a370462dfc5317cc6b207b31cc))

* docs: add sphinx itself to RTD requirements ([`44dbe7c`](https://github.com/lanl/ThunderBoltz/commit/44dbe7cdb0ee8d4a41d3d9fba16ea920d67bb282))

* docs: add RTD req.txt ([`4c8ce66`](https://github.com/lanl/ThunderBoltz/commit/4c8ce66f030ea28f0d17259f103d4b8179da2394))

* docs: remove RTD extension ([`ab41489`](https://github.com/lanl/ThunderBoltz/commit/ab41489aece7a71dbf1663d3cf49d9cc6b78d567))

* docs: remove rtd from requirements ([`cb83758`](https://github.com/lanl/ThunderBoltz/commit/cb83758fda716899f3f45c5938a73789ed7a56fc))

* docs: use ReadTheDocs theme ([`9996ffa`](https://github.com/lanl/ThunderBoltz/commit/9996ffa0ec774f939cdabdbe52b42bf3fa1d281c))

* docs: add requirements for ReadTheDocs ([`0c05d69`](https://github.com/lanl/ThunderBoltz/commit/0c05d6947165c3eeb14aac21cc305a054ff0b744))

* docs: add readthedocs workflow ([`c15fe16`](https://github.com/lanl/ThunderBoltz/commit/c15fe162844fd0e4acdd533cd58538e20c4c5ea2))

### Fix

* fix(docs): include custom css ([`70c45d1`](https://github.com/lanl/ThunderBoltz/commit/70c45d117a607ca81f231fb9bd06f23b6ff88dfd))


## v0.0.0 (2023-11-06)

### Build

* build: move workflow to its own dir ([`6d0759e`](https://github.com/lanl/ThunderBoltz/commit/6d0759eaa7017cb2ac2b30747bf386cc539252d9))

* build: rename action workflows ([`d92dcaf`](https://github.com/lanl/ThunderBoltz/commit/d92dcaf2f445065f0ab3a63319620fa1b66281f9))

* build: transition dev workflow to poetry / CICD ([`e0336ab`](https://github.com/lanl/ThunderBoltz/commit/e0336abab3889f80269a74e983c390ff5b07c614))

### Documentation

* docs: add thunderBoltz compilation clarity ([`42fe055`](https://github.com/lanl/ThunderBoltz/commit/42fe05554c767018b7d24e8e552033298806cf49))

* docs: Remove extraneous LA-UR ([`54eeb3e`](https://github.com/lanl/ThunderBoltz/commit/54eeb3e4c3a4572fc20cd679c8747822939e2944))

### Test

* test: add duration and verbose in pytest action ([`2dcb05d`](https://github.com/lanl/ThunderBoltz/commit/2dcb05dec3cdf8c19a0089d15d675cf993be64b8))

### Unknown

* Update README.md ([`76e8a0b`](https://github.com/lanl/ThunderBoltz/commit/76e8a0bc3843868575708294aa25d86f68cf9512))

* Add LA-URs to manuals and change titles ([`9bad8e1`](https://github.com/lanl/ThunderBoltz/commit/9bad8e1572788927b95e6717f11aa0c7809817c2))

* Remove CO collision ordering option ([`bdcb85e`](https://github.com/lanl/ThunderBoltz/commit/bdcb85e27e4dcf8556b30a11d38ae91d95877b2a))

* Update README with new cpp source location ([`4a2b826`](https://github.com/lanl/ThunderBoltz/commit/4a2b826cf8b32be0f911ef53c28fa49c3fa150fd))

* Ignore doc file outputs ([`40ba7d9`](https://github.com/lanl/ThunderBoltz/commit/40ba7d932a7cae4b2bf68236bba6d5d835e09f64))

* Rename pytb -&gt; thunderboltz and remove standalone ([`c2de216`](https://github.com/lanl/ThunderBoltz/commit/c2de216b0ea26217d1112577dae81e20a7032f68))

* ThunderBoltz v0.1 : Pre-paper

ThunderBoltz first commit for pre-paper submission.
To be cleaned (file/folder name changes, compiler friendly changes). ([`bdd3013`](https://github.com/lanl/ThunderBoltz/commit/bdd3013da1954440ed68eec30611d6dad479b6b3))

* Update README.md ([`70cef3a`](https://github.com/lanl/ThunderBoltz/commit/70cef3af613dc2cb4fc4ebf8547f61b39bf90790))

* Update README.md ([`2324cf4`](https://github.com/lanl/ThunderBoltz/commit/2324cf4fa11eaabb025c607dd4359e39970bfb9c))

* Update README.md ([`e4b0e02`](https://github.com/lanl/ThunderBoltz/commit/e4b0e02cd0ab82f471abd6b2fe8614f9de935971))

* Update README.md ([`cccac37`](https://github.com/lanl/ThunderBoltz/commit/cccac37c177396e0d4f7ad84d55b828a08d24bf4))

* Update README.md ([`324dfe9`](https://github.com/lanl/ThunderBoltz/commit/324dfe9cfcf2054b44bd2b5cda74caf90019b490))

* Update README.md ([`697b08f`](https://github.com/lanl/ThunderBoltz/commit/697b08f8af15d18b91c3b69c9a35a0d798e068b1))

* Update README.md ([`05745c4`](https://github.com/lanl/ThunderBoltz/commit/05745c4200f7392098ff5429bcba9c6c919219ec))

* Update README.md ([`17e515b`](https://github.com/lanl/ThunderBoltz/commit/17e515be10692cff58095841b6c8bb52bb1e488c))

* Update README.md ([`3fd1e69`](https://github.com/lanl/ThunderBoltz/commit/3fd1e694f666506cdec41a331af0f0696bd454fb))

* Update README.md ([`a875a2b`](https://github.com/lanl/ThunderBoltz/commit/a875a2bd2ecca6e066d0a376147ae093d5e0fbd1))

* Initial commit ([`322cbc7`](https://github.com/lanl/ThunderBoltz/commit/322cbc74d07c0eb69dc7fb400467003b3915f50d))
