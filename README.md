# Hyperpolarizability.com

This is the source code for [hyperpolarizability.com](https://hyperpolarizability.com), a website about quantum chemistry, computational physics, and related topics.

## Site Structure

This site is built with [Hakyll](https://jaspervdj.be/hakyll/), a static site generator written in Haskell.

- `/posts/`: Markdown files for blog posts
- `/templates/`: HTML templates
- `/css/`: CSS stylesheets
- `/js/`: JavaScript files
- `/images/`: Image files
- `/bib/`: BibTeX bibliography files for academic citations
- `site.hs`: Hakyll configuration file

## Local Development

### Prerequisites

- [GHC](https://www.haskell.org/ghc/) (Glasgow Haskell Compiler)
- [Cabal](https://www.haskell.org/cabal/) (Common Architecture for Building Applications and Libraries)

### Building the Site

1. Clone this repository:
   ```
   git clone https://github.com/noprofits-org/hyperpolarizability.com.git
   cd hyperpolarizability.com
   ```

2. Build the site generator:
   ```
   cabal build
   ```

3. Generate the site:
   ```
   cabal exec site build
   ```

4. Preview the site:
   ```
   cabal exec site watch
   ```
   
   This will start a local server at [http://localhost:8000](http://localhost:8000).

### Adding a New Post

1. Create a new Markdown file in the `/posts/` directory
2. Include YAML frontmatter with title, date, and other metadata
3. Write your content using Markdown
4. For academic citations, use the format `[@citation-key]` to reference entries in the bibliography

Example post frontmatter:
```markdown
---
title: My New Post
date: 2025-04-19
tags: tag1, tag2
description: Brief description of the post
bibliography: /path/to/bibliography.bib
csl: /path/to/style.csl
---

Post content goes here...
```

## Deployment

The site is automatically deployed to GitHub Pages when changes are pushed to the `main` branch using GitHub Actions. The deployment workflow:

1. Builds the Hakyll site generator
2. Generates the static site
3. Deploys the generated files to the `gh-pages` branch
4. GitHub Pages serves the content from this branch

## License

See the [LICENSE](LICENSE) file for details.