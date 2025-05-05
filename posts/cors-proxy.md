---
title: Building a Serverless CORS Proxy with Vercel: Simplifying Cross-Origin Requests
date: 2025-05-04
tags: cors proxy, serverless, vercel, javascript, web development
description: A practical exploration of building a serverless CORS proxy using Vercel's serverless functions, offering an elegant solution to the common cross-origin resource sharing challenges faced by frontend developers.
---

## Cross-Origin Headaches and Their Serverless Solution

Cross-origin resource sharing (CORS) remains one of those persistent challenges that frontend developers regularly encounter. You build a beautiful client-side application, everything works perfectly in development, and then you deploy it only to see those dreaded red error messages: "Access to fetch at 'https://api.example.com' from origin 'https://your-app.com' has been blocked by CORS policy." This frustration is especially common when attempting to access third-party APIs that haven't configured CORS to allow requests from your domain.

After encountering this issue repeatedly in our nonprofit applications, I decided to create a simple yet effective solution: a serverless CORS proxy deployed on Vercel. Our proxy server acts as a middleman between your client-side code and any API you want to access, effectively bypassing CORS restrictions by making server-side requests on your behalf.

## How Our CORS Proxy Works

The concept is straightforward: instead of making a direct request from your browser to a target API (which would be blocked by CORS restrictions), you send your request to our proxy server. The proxy then forwards your request to the target API, receives the response, and returns it to your application with the appropriate CORS headers that allow your browser to accept the response.

The implementation leverages Vercel's serverless functions, which provide a lightweight and scalable infrastructure without the need to maintain dedicated servers. The core of our solution is a single JavaScript file (proxy.js) that handles incoming requests, forwards them to the specified target URL, and returns the response with properly configured CORS headers.

The most challenging aspect of the implementation was handling different content types correctly. APIs might return JSON data, HTML content, images, or other binary data. Our proxy needed to properly identify the content type and pass it through without corruption. Additionally, I needed to ensure that the proxy would work with various HTTP methods (GET, POST, PUT, DELETE, etc.) and correctly forward headers and request bodies.

After several iterations and testing, I developed a solution that handles these requirements elegantly. One key insight was to exclude certain headers that could cause conflicts or security issues (like host, connection, origin, and content-encoding) when forwarding requests.

## Using the CORS Proxy in Your Projects

Integrating our CORS proxy into your projects is remarkably simple. For a basic GET request, you would use:

```javascript
const PROXY_URL = 'https://cors-proxy-xi-ten.vercel.app/api/proxy';
const TARGET_API = 'https://api.example.com/data';

fetch(`${PROXY_URL}?url=${encodeURIComponent(TARGET_API)}`)
  .then(response => response.json())
  .then(data => {
    console.log('Data:', data);
  })
  .catch(error => console.error('Error:', error));
```

For POST requests or other methods that require a body and specific headers:

```javascript
const PROXY_URL = 'https://cors-proxy-xi-ten.vercel.app/api/proxy';
const TARGET_API = 'https://api.example.com/data';

fetch(`${PROXY_URL}?url=${encodeURIComponent(TARGET_API)}`, {
  method: 'POST',
  headers: {
    'Content-Type': 'application/json'
  },
  body: JSON.stringify({ key: 'value' })
})
  .then(response => response.json())
  .then(data => {
    console.log('Data:', data);
  })
  .catch(error => console.error('Error:', error));
```

The proxy seamlessly handles different content types, including JSON, text, and binary data like images, making it versatile for various API integrations.

## Testing the Proxy

To help developers test the proxy with different request configurations, I created a dedicated testing page at [https://noprofits.org/cors-tester/](https://noprofits.org/cors-tester/). This interactive tool allows you to:

1. Enter a target URL to access through the proxy
2. Select from various HTTP methods (GET, POST, PUT, etc.)
3. Configure request headers and body content
4. Execute the request and view the response in a user-friendly format

The tester provides real-time feedback and displays the response data, headers, and any errors that might occur. It's particularly useful for debugging API interactions and understanding how different request configurations affect the results.

## Technical Challenges and Solutions

Building the proxy wasn't without challenges. One particularly troublesome issue was handling content encoding properly. I initially encountered "ERR_CONTENT_DECODING_FAILED" errors when proxying certain responses. After investigation, I discovered that this was happening because the content-encoding header was being forwarded from the target API to the client, but the content had already been decoded by the node-fetch library I were using.

The solution was to explicitly remove the content-encoding header from responses before sending them back to the client:

```javascript
// Explicitly remove content-encoding header to prevent decoding errors
res.removeHeader('content-encoding');
```

Another challenge was handling various content types correctly. For JSON responses, I needed to parse and re-stringify the content. For binary data like images, I needed to forward the raw buffer without modification:

```javascript
// Handle different content types appropriately
const contentType = response.headers.get('content-type') || '';

if (contentType.includes('application/json')) {
  const jsonData = await response.json();
  return res.json(jsonData);
} else if (contentType.includes('text/')) {
  const text = await response.text();
  return res.send(text);
} else {
  // Handle binary data (e.g., images)
  const buffer = await response.buffer();
  return res.send(buffer);
}
```

## Open Source and Available for Everyone

We've made our CORS proxy open source and available on GitHub at [https://github.com/noprofits-org/cors-proxy-server](https://github.com/noprofits-org/cors-proxy-server). Developers can either use our hosted version directly or fork the repository to deploy their own instance on Vercel.

The repository includes comprehensive documentation, example code, and deployment instructions. It's designed to be easy to understand and modify, even for developers who aren't familiar with serverless functions or CORS concepts.

## Beyond Development: Production Considerations

While our CORS proxy is a valuable tool for development and testing, there are important considerations for production use. The proxy is currently open without authentication, making it accessible to anyone. For production applications with significant traffic or sensitive data, I recommend deploying your own instance with additional security measures.

The current implementation relies on Vercel's generous free tier, which includes reasonable limits for most small to medium-sized applications. However, very high-traffic applications might require upgrading to a paid plan or implementing more sophisticated caching and rate-limiting strategies.

## Conclusion

Cross-origin resource sharing doesn't have to be a roadblock in your development process. With our serverless CORS proxy, you can easily bypass these restrictions and focus on building great applications rather than wrestling with browser security policies.

I built this tool to solve our own challenges when developing applications for nonprofits, and we're excited to share it with the broader development community. Whether you're building a small personal project or a complex application, our CORS proxy provides a simple, effective solution to a common problem.

Visit [https://noprofits.org/cors-tester/](https://noprofits.org/cors-tester/) to try the proxy for yourself, and check out the [GitHub repository](https://github.com/noprofits-org/cors-proxy-server) to learn more about the implementation or deploy your own instance. Happy coding!