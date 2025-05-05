/**
 * print.js - Handles print functionality
 * 
 * This script adds print functionality to the print button.
 * It provides a clean print experience for the user.
 */

document.addEventListener('DOMContentLoaded', function() {
  // Find all print buttons
  const printButtons = document.querySelectorAll('.print-button');
  
  // Add click event listener to each print button
  printButtons.forEach(button => {
    button.addEventListener('click', function(e) {
      e.preventDefault();
      
      // Trigger the print dialog
      window.print();
    });
  });
});