/**
 * code-copy.js - Adds copy buttons to code blocks
 * 
 * This script adds copy-to-clipboard functionality to all code blocks.
 * It creates a button that appears when hovering over code blocks,
 * and copies the code content when clicked.
 */

document.addEventListener('DOMContentLoaded', function() {
  // Find all pre elements (code blocks)
  const codeBlocks = document.querySelectorAll('pre');
  
  // For each code block, add a copy button
  codeBlocks.forEach((block) => {
    // Create button element
    const button = document.createElement('button');
    button.className = 'copy-button';
    button.textContent = 'Copy';
    
    // Add the button to the code block
    block.appendChild(button);
    
    // Add click event listener to the button
    button.addEventListener('click', () => {
      // Get the text content of the code element
      const code = block.querySelector('code') 
        ? block.querySelector('code').textContent
        : block.textContent;
      
      // Copy to clipboard
      navigator.clipboard.writeText(code)
        .then(() => {
          // Visual feedback for success
          button.textContent = 'Copied!';
          button.style.backgroundColor = 'var(--accent-primary)';
          button.style.color = 'white';
          
          // Reset after 2 seconds
          setTimeout(() => {
            button.textContent = 'Copy';
            button.style.backgroundColor = '';
            button.style.color = '';
          }, 2000);
        })
        .catch(err => {
          // Visual feedback for failure
          console.error('Failed to copy: ', err);
          button.textContent = 'Failed';
          button.style.backgroundColor = '#f44336';
          button.style.color = 'white';
          
          // Reset after 2 seconds
          setTimeout(() => {
            button.textContent = 'Copy';
            button.style.backgroundColor = '';
            button.style.color = '';
          }, 2000);
        });
    });
  });
});